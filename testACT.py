from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
import healpy as hp
import numpy as np

ACT_148_table = fits.open('../catalogs/ACT_148.fits')
data = ACT_148_table[0].data
head = ACT_148_table[0].header

RA = []
DEC = []
for i in range(head['NAXIS1']):
    ra_val = (i-head['CRPIX1'])*head['CDELT1']
    RA.append(ra_val)
    
for i in range(head['NAXIS2']):
    dec_val = (i-head['CRPIX2'])*head['CDELT2']
    DEC.append(dec_val)
    
#Make NxM arrays of ra and dec this way!
data_ra = np.zeros([head['NAXIS2'], head['NAXIS1']])
data_dec = np.zeros([head['NAXIS2'], head['NAXIS1']])

print('beginning calculating ra / dec')

for i in range(head['NAXIS2']):
    for j in range(head['NAXIS1']):
        dec_val = (i-head['CRPIX2'])*head['CDELT2']
        ra_val = (j-head['CRPIX1'])*head['CDELT1']
        data_ra[i, j] = ra_val
        data_dec[i, j] = dec_val

print('finished with calculating ra / dec')

hdu_data = fits.PrimaryHDU(data)
hdu_ra = fits.ImageHDU(data_ra)
hdu_dec = fits.ImageHDU(data_dec)

hdu_list = fits.HDUList([hdu_data, hdu_ra, hdu_dec])
hdu_list.writeto('../catalogs/ACT_148_withCoord.fits')

print('Min corner:')
print(data_ra[0, 0])
print(data_dec[0,0])
print('Center:')
print(data_ra[292, 7972])
print(data_dec[292, 7972])
print('Max corner:')
print(data_ra[-1, -1])
print(data_dec[-1, -1])
