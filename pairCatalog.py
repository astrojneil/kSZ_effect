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
#ind_i, ind_j, r_ij, theta_ij, cij
column_list=['ind_i', 'ind_j', 'r_ij', 'theta_ij', 'c_ij']
df = pd.DataFrame(columns=column_list)


print('Beginning pair calculation....')
d, p = initBar()
#scales as (len(temp_data)*len(temp_data)-1)/2
for i in r: #range(len(temp_data)):
    c1 = SkyCoord(ra = temp_data['ra'][i]*u.degree, dec =temp_data['dec'][i]*u.degree, distance = temp_data['R'][i]*u.Mpc)
    for j in r:  #range(len(temp_data)):
        if i < j:
            c2 = SkyCoord(ra = temp_data['ra'][j]*u.degree, dec = temp_data['dec'][j]*u.degree, distance= temp_data['R'][j]*u.Mpc)
            theta = c1.separation(c2).rad
            dist = c1.separation_3d(c2)
            r_i = temp_data['R'][i]
            r_j = temp_data['R'][j]

            #skip over pairs with separations greater than 400 Mpc in order to reduce final pair catalog size
            if dist < 400*u.Mpc:
                c_ij = ((r_i - r_j)*(1+math.cos(theta)))/(2*np.sqrt(r_i**2 + r_j**2 - 2*r_i*r_j*math.cos(theta)))

                newrow = [i, j, dist.value, theta, c_ij]
                newrow_frame = pd.DataFrame([newrow], columns = column_list)
                df = pd.concat([df,newrow_frame], ignore_index = True)
            else:
                #skip this pair
                continue
        else:
            #avoid over-counting
            continue
    d = d+1
    d, p = updateBar(d, p, len(temp_data))

print('\nDone!')
df.to_csv('pairCatalog_ind_short.csv')
