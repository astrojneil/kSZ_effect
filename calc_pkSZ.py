#This script calculates the pkSZ signal within each separation bin using the
#pairs of galaxies from the catalog created with pairCatalog.py
#The errors are estimated with a bootstrap resampling technique; resampling
#each distance bin 1000 times. The standard deviation of the distribution of 
#these samples is then quoted as the error on each point. 

#The results are then plotted pkSZ vs. Separation

#imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from progressBar import *
import random

#load data
pairData = pd.read_csv('pairCatalog.csv')

#define bins
binSize = 30
numBins = 10
binStart = 0

#define arrays to hold data
all_sums = []
pkSZ_all = []
binCenters = []
binCount = []
errs = []

#loop through bins
for i in range(numBins):
    binEnd = binStart+binSize
    binCen = binStart + binSize/2
    binCenters.append(binCen)

    print('Bin: {} to {}'.format(binStart, binEnd))
    
    #select rows for each bin
    bin_data = pairData[(pairData['r_ij'] >= binStart) & (pairData['r_ij'] < binEnd)]
    N = len(bin_data['r_ij'])
    binCount.append(N)

    #sum or 'stack' all pairs within the bin
    bot = bin_data['c_ij']*bin_data['c_ij']
    bot_sum = bot.sum()

    top = ((bin_data['Temp_1']-bin_data['T_z_1'])-(bin_data['Temp_2']-bin_data['T_z_2']))*bin_data['c_ij']
    top_sum = top.sum()

    pkSZ= -top_sum/bot_sum
    pkSZ_all.append(pkSZ)

    #calculate errors by resampling bin 1000 times
    binSums = []
    d, p = initBar()
    for j in range(1000):
        r = [random.randint(0, N-1) for iter in range(N)]
        newBin = bin_data.iloc[r]

        #sum or 'stack' all pairs within the bin
        bot = newBin['c_ij']*newBin['c_ij']
        bot_sum = bot.sum()

        top = ((newBin['Temp_1']-newBin['T_z_1'])-(newBin['Temp_2']-newBin['T_z_2']))*newBin['c_ij']
        top_sum = top.sum()

        pkSZ= -top_sum/bot_sum
        binSums.append(pkSZ)
        d = d+1
        d, p = updateBar(d, p, 1000)

    errs.append(np.std(binSums))  #save 1 sigma errors
    all_sums.append(binSums)
    binStart = binEnd
    print('\n- - - - - -  bin done  - - - - - -')

#plot results
fig = plt.figure()
ax = fig.add_subplot()

#pkSZ vs Sep
ax.errorbar(binCenters, pkSZ_all, yerr = errs, fmt = 'ro')
ax.plot(binCenters, np.zeros(len(binCenters)), color = 'black')
ax.set_ylabel('pkSZ (microKelvin)')
ax.set_xlabel('Separation (Mpc)')


#historgram:
#ax.hist(all_sums[0], color = 'red', bins = 50)
#ax.plot([pkSZ_all[0], pkSZ_all[0]], [0, 90], linestyle = 'dashed', color = 'black')
#ax.plot([pkSZ_all[0]+errs[0], pkSZ_all[0]+errs[0]], [0, 90], linestyle = 'dashed', color = 'black')
#ax.plot([pkSZ_all[0]-errs[0], pkSZ_all[0]-errs[0]], [0, 90], linestyle = 'dashed', color = 'black')

#save figure
fig.savefig('pkSZ_v_Sep.png')
