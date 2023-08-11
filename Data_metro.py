from astropy.io import fits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import random

"""
#MAYBE SOME OTHER TIME
import cosmoglobe 
model = cosmoglobe.sky_model(nside=256)
"""  

nside = 512

def signal_to_noise(N, band):
    
    filename = f'Galactic_2021-10-19_Feeds1-9-10-11-14-15-18-19_Band{band}.fits'
    savename = filename[:-5] #removes the .fits for the filename to be saved in a new format

    signal_map = hp.read_map(f'./COMAP_to_healpix/Fits_files/Signal_healpix_{savename}.fits')
    uncertainty_map = hp.read_map(f'./COMAP_to_healpix/Fits_files/Uncertainty_healpix_{savename}.fits')


    """Finding signal to noise ratio"""

    signal_map_ma = hp.ma(signal_map)
    uncertainty_map_ma = hp.ma(uncertainty_map)

    signal_to_noise = np.abs(signal_map_ma/uncertainty_map_ma)

    return signal_to_noise, signal_map_ma, uncertainty_map_ma

idx = 23646863

def sort_data(index):
    
    data = np.array([signal_to_noise(nside, i)[1][idx] for i in range(8)])
    uncertainty = np.array([signal_to_noise(nside, i)[2][idx] for i in range(8)])
    
    np.savetxt(f'metro_data_index_{idx}.txt', data, delimiter=',')
    np.savetxt(f'metro_uncertainty_index_{idx}.txt', uncertainty, delimiter=',')

    return

save_file = sort_data(idx)


