# import statements
import emcee
import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
import corner
import scipy.stats
import h5py
import csv
from astropy.io          import ascii


nburn = 80000

#read MCMC data for each snr
#snr 20:
snr20_reader = emcee.backends.HDFBackend('snr20.h5')
# get samples, discarding burn-in and apply thinning
snr20_samples  = snr20_reader.get_chain(discard = nburn)
snr20_ndim     = snr20_samples.shape[2]

#flatten the chain
flatchain_20 = snr20_samples.reshape((-1,snr20_ndim))

#land posteriors are 7 and 8

basalt_20 = 10**(flatchain_20[:,7]).reshape(-1,1)
granite_20 = 10**(flatchain_20[:,8]).reshape(-1,1)

sum_20 = (basalt_20 + granite_20)

#snr30
snr30_reader = emcee.backends.HDFBackend('snr30.h5')
# get samples, discarding burn-in and apply thinning
snr30_samples  = snr30_reader.get_chain(discard = nburn)
snr30_ndim     = snr30_samples.shape[2]

#flatten the chain
flatchain_30 = snr30_samples.reshape((-1,snr30_ndim))

basalt_30 = 10**(flatchain_30[:,7]).reshape(-1,1)
granite_30 = 10**(flatchain_30[:,8]).reshape(-1,1)

sum_30 = (basalt_30 + granite_30)

#snr40
snr40_reader  = emcee.backends.HDFBackend('snr40.h5')
snr40_samples = snr40_reader.get_chain(discard = nburn)
snr40_ndim    = snr40_samples.shape[2]

#flatten chain
flatchain_40 = snr40_samples.reshape((-1, snr40_ndim))

basalt_40 = 10**(flatchain_40[:,7]).reshape(-1,1)
granite_40 = 10**(flatchain_40[:,8]).reshape(-1,1)

sum_40 = (basalt_40 + granite_40).reshape(-1,1)

#make a chain with the individual land fractions and the total

full_40 = np.stack((basalt_40[:,0], granite_40[:,0], sum_40[:,0]), axis = 1)
full_30 = np.stack((basalt_30[:,0], granite_30[:,0], sum_30[:,0]), axis = 1)
full_20 = np.stack((basalt_20[:,0], granite_20[:,0], sum_20[:,0]), axis = 1)

figure = corner.corner(full_40, color = 'purple', labels = ['W basalt', 'Granite', 'Total land'], truths = [0.1,0.2, 0.3], truth_color = 'black', plot_contours = True)
corner.corner(full_30, color = 'green', labels = ['W basalt', 'Granite', 'Total land'], truths = [0.1,0.2, 0.3], truth_color = 'black' , plot_contours = True, fig = figure)   
corner.corner(full_20, color = 'blue', labels = ['W basalt', 'Granite', 'Total land'], truths = [0.1,0.2,0.3], truth_color = 'black', plot_contours = True, fig = figure)    
plt.figlegend(['SNR 40', 'Truth = 0.1', 'SNR 30', 'Truth = 0.2', 'SNR 20', 'Truth = 0.3'])
figure.savefig('land_corner.png',format='png',bbox_inches='tight', dpi = 300)
plt.close()





