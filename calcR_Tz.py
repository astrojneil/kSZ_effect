#this script calculates two new columns for the boss_cleaned
#dataframe:
#r, the comoving distance of the galaxy
#t_z = the contribution to the temperature of a single galaxy, i, from all
#other galaxies, with the galaxies at redshifts near z_i weighted the most strongly
#according to a gaussian centered on z_i
#These columns are then saved the boss_cleaned dataframe, and saved

#imports
import numpy as np
import pandas as pd
from progressBar import *

#constants, Hand et al., Dunkle et al. 2011
H_0 = 69.7  #km/s/Mpc
om_m = 0.279
c = 3e5   #km/s
sig_z = 0.01
R_cons = (2*c)/(H_0*np.sqrt(om_m))

#open file
boss_clean = pd.read_csv('cleaned_boss_temps.csv')

RList = []
T_zList = []
print('Beginning calculations...')
i = 0
j = 0
initBar()
for g in range(len(boss_clean)):
    R = R_cons*(1 - 1/(np.sqrt(1+boss_clean['z'][g])))
    RList.append(R)

    T_bot = np.exp(-(boss_clean['z'][g] - boss_clean['z'])**2/(2*sig_z**2))
    T_bot_sum = T_bot.sum()

    T_top = boss_clean['CMB_temp']*np.exp(-(boss_clean['z'][g] - boss_clean['z'])**2/(2*sig_z**2))
    T_top_sum = T_top.sum()
    
    T_z = T_top_sum/T_bot_sum
    T_zList.append(T_z)

    i = i+1
    i, j = updateBar(i, j, len(boss_clean))


print('\nDone! Saving to catalog')
boss_clean['R'] = RList
boss_clean['T_z'] = T_zList

boss_clean.to_csv('cleaned_boss_Temp_R_Tz.csv')