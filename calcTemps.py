#This script takes ACT data and a cleaned BOSS catalog and computes
#the average temperature within 1 arcminute centered on a galaxy within BOSS
#ACT data is sampled as 10 arcminute boxes and repixelized to 0.065'
#The repixelized map is convolved with the ACT beam and normalized before
#masked and averaged over the 1 arcminute ceneterd on the galaxy.
#Calculated temperatures are saved as a column in the cleaned_boss catalog

#imports
from astropy.io import fits
from astropy import wcs
import pandas as pd
import numpy as np
from scipy import ndimage
from scipy import interpolate
from progressBar import *
import numpy.ma as ma
from scipy.integrate import simps


#function to mask data within a radius of a certain point
def datamask(data, radius, point):
    x = np.linspace(0,len(data)-1, len(data))
    y = np.linspace(0, len(data[0])-1, len(data[0]))

    xx, yy = np.meshgrid(x, y)
    zz = np.zeros((len(data), len(data[0])))
    size = 0
    for i in range(len(x)):
        for j in range(len(y)):
            x_ind = xx[i][j]+1
            y_ind = yy[i][j]+1
            rad = np.sqrt((x_ind-point[0])**2+(y_ind-point[1])**2)
            if rad < radius:
                zz[i][j] = 1
                size = size +1

    return zz, size

#open files
print('Opening Files and Initializing...')
ACT_table = fits.open('../../catalogs/ACT_148.fits')
ACT_hits = fits.open('../../catalogs/ACT_148_equ_season_3_1way_hits_v3.fits')
ACT_data = ACT_table[0].data
ACT_header = ACT_table[0].header
hits = ACT_hits[0].data

boss_clean = pd.read_csv('../../catalogs/cleaned_boss.csv')

#define pixel/world coordinates
w = wcs.WCS(ACT_header)

#constants for map making and repixelizing
init_res_x = 20
init_res_y = 20
box_x = init_res_x-1
box_y = init_res_y-1
new_res = 167  #0.0625 is 1/16 of an arcminute... so there are 160 cells in 10 arcminutes
rad = 16  # 1 arcminute in pixels at the new resolution

x = np.linspace(0, box_x, init_res_x)
y = np.linspace(0, box_y, init_res_y)
xx, yy = np.meshgrid(x, y)

xnew = np.linspace(0, box_x, new_res)
ynew = np.linspace(0, box_y, new_res)

#define ACT beam as a 2D array by mirroring the profile across 0 degrees
#only need to include the radii out to 0.0083 (5 arcminutes) (first 83 items)
#this also ensures that the beam profile is the same dimensions as the
#newly pixelized map
radius = []
beam = []

beamfile = open('../../catalogs/profile_AR1_2008_pixwin_130224.txt')

#read in beam data
for line in beamfile:
    ls = line.split()
    radius.append(float(ls[0]))
    beam.append(float(ls[1]))

radius = np.array(radius)
beam = np.array(beam)
beamfile.close()

flipBeam = np.zeros(83*2+1)
extraRad = np.zeros(83*2+1)

for i in range(83*2+1):
    if (i-83) <= 0:
        flipBeam[i] = beam[83-i]
        extraRad[i] = -0.0083+radius[i]
    else:
        extraRad[i] = radius[i-83]
        flipBeam[i] = beam[i-83]

#multipy this beam profile by itself as a vector to create a 2D profile of the beam
flipMat = np.matrix(flipBeam)
beam2D = np.asarray(np.multiply(flipMat, flipMat.T))

#normalization for beam convolution
cert_beam = np.ones((new_res, new_res))
norm_beam = ndimage.convolve(cert_beam, beam2D, mode='constant', cval = 1.0)

#find the pixel at the location of each BOSS galaxy
list_coord = []
for i in range(len(boss_clean)):
    ra = boss_clean['ra'][i]
    dec = boss_clean['dec'][i]

    coord = [ra, dec]
    list_coord.append(coord)

center_pix = w.wcs_world2pix(list_coord, 1)
print(np.average(ACT_data))

#create matched filter;
weights = np.sqrt(hits/np.max(hits))
newdata = weights*ACT_data
print(np.average(newdata))

#T_other = data with boss sources masked
#make the mask
#also set certainty array to zero where sources are masked
#(there is no signal to the filter from masked locations)
#(but does the gauassian smooth part put signal there?)
mask = np.zeros((len(ACT_data), len(ACT_data[0])))
cert_map = np.ones((len(ACT_data), len(ACT_data[0])))
#loop through sources
for c in range(len(center_pix)):
    dec_pix = int(round(center_pix[c,1]))
    ra_pix = int(round(center_pix[c,0]))

    mask[dec_pix-3:dec_pix+3, ra_pix-3:ra_pix+3] = 1
    cert_map [dec_pix-3:dec_pix+3, ra_pix-3:ra_pix+3] = 0


T_other = ma.masked_array(ACT_data, mask=mask)
T_other_FT = np.fft.fft2(T_other)

T_other_abs = np.absolute(T_other_FT)**(-2)
#smooth T_other with a gaussian
T_other_FT_smooth = ndimage.gaussian_filter(T_other_abs, 5)

#find median value
med = np.median(T_other_FT_smooth)
#set any pixel > 10*median to the median value
T_other_FT_smooth[T_other_FT_smooth > 10*med]

#beam FT
beam_FT_2d = np.fft.fft2(beam2D)
beam_conj = np.conj(beam_FT_2d)

norm_beamFFT = ndimage.convolve(cert_map, beam_FT_2d, mode='constant', cval=1.0)
norm_beam_conj = ndimage.convolve(cert_map, beam_conj, mode='constant', cval=1.0)

#filter
filt_top = ndimage.convolve(T_other_FT_smooth, beam_conj, mode= 'constant', cval = 0.0)
filt_top = filt_top/norm_beam_conj
bot1 = ndimage.convolve(filt_top, beam_FT_2d, mode= 'constant', cval = 0.0)
bot1 = bot1/norm_beamFFT

x_int = np.linspace(0, 1, len(bot1))
y_int = np.linspace(0, 1, len(bot1[0]))
#integrate filt_bottom to get the normalizing factor
filt_bottom = simps(simps(bot1, y_int), x_int)
match_filt = filt_top/filt_bottom

#apply the filter to the weighted data
finaldata = np.abs(np.fft.ifft2(np.fft.fft2(newdata)*match_filt))
print(np.average(finaldata))

#loop through BOSS galaxies, pull out 10' map around the center pixel
#repixelize, and convolve the map, saving the average temperatures
print('Beginning loop through BOSS galaxies:')

all_temps = []

#progress bar:
i,j = initBar()

for c in range(len(boss_clean)):
    ra_pix = int(round(center_pix[c,0]))
    dec_pix = int(round(center_pix[c,1]))

    dataMap = finaldata[dec_pix-10:dec_pix+10, ra_pix-10:ra_pix+10]
    interp = interpolate.interp2d(x, y, dataMap)
    newMap = interp(xnew, ynew)

    convMap = ndimage.convolve(newMap, beam2D, mode='constant', cval = 0.0)
    normMap = convMap/norm_beam

    maskAr, size = datamask(normMap, 16, [int(new_res/2), int(new_res/2)])
    temp = np.sum(maskAr*normMap)/size
    all_temps.append(temp)
    i = i+1
    i, j = updateBar(i, j, len(boss_clean))



print('Saving temperatures to catalog')
boss_clean['CMB_temp'] = all_temps

boss_clean.to_csv('cleaned_boss_temps_test.csv')

ACT_table.close()
