import numpy as np
import pandas as pd
from progressBar import *
import matplotlib.pyplot as plt


pairData = pd.read_csv('pairCatalog.csv')


#define bins
binSize = 9
numBins = 20
binStart = 0


pkSZ_all = []
binCenters = []
binCount = []
binErr = []
for i in range(numBins):
    binEnd = binStart+binSize
    binCen = binStart + binSize/2
    binCenters.append(binCen)

    #select rows for each bin
    bin_data = pairData[(pairData['r_ij'] >= binStart) & (pairData['r_ij'] < binEnd)]
    binCount.append(len(bin_data['r_ij']))
    #sum or 'stack' all pairs within the bin
    bot = bin_data['c_ij']*bin_data['c_ij']
    bot_sum = bot.sum()

    top = ((bin_data['Temp_1']-bin_data['T_z_1'])-(bin_data['Temp_2']-bin_data['T_z_2']))*bin_data['c_ij']
    top_sum = top.sum()

    pkSZ= top_sum/bot_sum
    pkSZ_all.append(pkSZ)
    binStart = binEnd

lastBin = pairData[(pairData['r_ij'] >= binStart)]
print("Unused pairs: {}".format(len(lastBin['r_ij'])))
print("Bin sizes:")
print(binCount)

print("Total Used Pairs: {} ({:.1%})".format(np.sum(np.array(binCount)), np.sum(np.array(binCount))/len(lastBin['r_ij']) ))

fig = plt.figure()
ax = fig.add_subplot()
ax.scatter(binCenters, pkSZ_all, marker = 'o', color = 'red')
ax.plot(binCenters, np.zeros(len(binCenters)), color = 'black')

fig.savefig('test_pkSZ.png')
