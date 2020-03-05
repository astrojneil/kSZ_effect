#This script creates a new catalog with the CMB temperatures, separation and distances
#for pairs of galaxies within the BOSS catalog
#This is done by looping through all galaxies, i, and calculating distance
#between all other galaxies, j, such that i < j to avoid over-counting.

#imports
import numpy as np
import pandas as pd
import math
from astropy.coordinates import SkyCoord
from astropy import units as u
from progressBar import *

import random

#load data
all_data = pd.read_csv('cleaned_boss_Temp_R_Tz.csv')

#randomly select 300 galaxies to consider average quantities
r = random.sample(list(np.arange(len(all_data))), k=300)
temp_data = all_data.iloc[r]

#create new catalog for pairs
#galaxy i: ra, dec, z, Temp, T_z, R
#galaxy j: ra, dec, z, Temp, T_z, R
#r_ij, theta_ij, pkSZ
column_list=['ra_1', 'dec_1', 'z_1', 'Temp_1', 'T_z_1', 'R_1','ra_2', 'dec_2', 'z_2', 'Temp_2', 'T_z_2', 'R_2', 'r_ij', 'theta_ij', 'c_ij']
df = pd.DataFrame(columns=column_list)


print('Beginning pair calculation....')
d, p = initBar()
#scales as (len(temp_data)*len(temp_data)-1)/2
for i in r: #range(len(temp_data)):
    c1 = SkyCoord(ra = temp_data['ra'][i]*u.degree, dec =temp_data['dec'][i]*u.degree)
    for j in r: #range(len(temp_data)):
        if i < j:
            r_i = temp_data['R'][i]
            r_j = temp_data['R'][j]
            r_ij = abs(r_i - r_j)

            #skip over pairs with separations greater than 750 Mpc in order to reduce final pair catalog size
            c2 = SkyCoord(ra = temp_data['ra'][j]*u.degree, dec = temp_data['dec'][j]*u.degree)
            theta = c1.separation(c2).rad  #math.cos expects radians

            c_ij = ((r_i - r_j)*(1+math.cos(theta)))/(2*np.sqrt(r_i**2 + r_j**2 - 2*r_i*r_j*math.cos(theta)))

            newrow = [temp_data['ra'][i], temp_data['dec'][i], temp_data['z'][i], temp_data['CMB_temp'][i], temp_data['T_z'][i], temp_data['R'][i],
                        temp_data['ra'][j], temp_data['dec'][j], temp_data['z'][j], temp_data['CMB_temp'][j], temp_data['T_z'][j], temp_data['R'][j],
                        r_ij, theta, c_ij]
            newrow_frame = pd.DataFrame([newrow], columns = column_list)
            df = pd.concat([df,newrow_frame], ignore_index = True)
        else:
            #avoid over-counting
            continue
    d = d+1
    d, p = updateBar(d, p, len(temp_data))

print('\nDone!')
df.to_csv('pairCatalog.csv')
