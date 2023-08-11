from astropy.io import fits
import numpy as np
import healpy as hp
import reproject  
import matplotlib.pyplot as plt
import sys
import time

"""
Unfinished
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

    #index = np.where(signal_to_noise>10)[0]  #Removing sn under 10
    #print(f'THE INDEX YOU ARE LOOKING FOR IS {index}')
    # signal_map_ = hp.ma(signal_map_ma)
    # #print(signal_map_)
    # uncertainty_map_ma = hp.ma(uncertainty_map_ma)
    # signal_to_noise = hp.ma(signal_to_noise)

    return signal_to_noise, signal_map_ma, uncertainty_map_ma




"""Picking a pixel which will remain the same always"""

# PICKING BAND 3
tol = 10
center = 35
mini = 30

sn_pix = (signal_to_noise(nside, 1))[0]

pixels = [np.intersect1d(np.where(sn_pix > mini)[0], np.where(np.isclose(sn_pix, center, atol=tol))[0])] #cutting 1 and assuming this works

ooga_boogas = pixels[0] #np.intersect1d(pixels[0], pixels[1]) #np.intersect1d(ooga, boogas)

#print(ooga_boogas[0])



"""Run through all bands and plot pixel brightness over the 8 bands"""
#Running through all the bands       ##ASSUMING IF ONE PIXEL IN ONE BAND HAS A SN>30, THE SAME PIXEL IN OTHER BANDS ALSO HAS SN>30##

# TIME 102s 
start = time.time()
sn = [signal_to_noise(nside, i)[1] for i in range(8)]
un = [signal_to_noise(nside, i)[2] for i in range(8)]

m_list = [e[ooga_boogas] for e in sn]
un_m_list = [u[ooga_boogas] for u in un]
# plt.show()
end = time.time()

for i in range(8):
    plt.errorbar(i, m_list[i][0], yerr = un_m_list[i][0], fmt = '*')  #picked the first
plt.title(f'Random pixel selcted with values of min {mini}, centered at {center} with a range of {tol}')
plt.xlabel('Band')
plt.ylabel('Signal')
plt.show()



""" Spectral index plotting ***"""

# TIME 98 s

m_list = []
sn_map_list = []

# start = time.time()

for i in range(8):
    sn = signal_to_noise(nside, i)[1]
    un = signal_to_noise(nside, i)[2]
    m_list.append(sn[ooga_boogas][0])
    sn_map_list.append(sn)
    
    plt.scatter(i, sn[ooga_boogas][0])  #picked the first
    plt.errorbar(i, sn[ooga_boogas][0], yerr = un[ooga_boogas][0], fmt = '*')  #picked the first
plt.title(f'Random pixel selcted with values of min {mini}, centered at {center} with a range of {tol}')
plt.xlabel('Band')
plt.ylabel('Signal')
plt.show()

# end = time.time()


"""Finding spectral index"""

m1 = m_list[-2]
m2 = m_list[-1]

nu1 = 32.5 #GHz  #freq for band[-2]
nu2 = 33.5 #GHz  #freq for band[-1]

beta = (np.log(m2)-np.log(m1))/(np.log(nu2)-np.log(nu1))
beta1 = np.log(m1)/np.log(nu1)
beta2 = np.log(m2)/np.log(nu2)

# print(beta)  

# print(beta1)

# print(beta2)
"""Output beta: -2.26096093876478"""

"""Make beta map, excluding SN>10""" 

#NOW MAKING BETA MAP

def maps_masked(band1, band2, noise, nside):   #Returns betamap, and signal maps used in the beta mapping

    nu = [26.7, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5]

    signal_map1 =  (signal_to_noise(nside, band1))[1]
    m1_map_snr = (signal_to_noise(nside, band1))[0]
    m1 = signal_map1.copy()
    m1[m1_map_snr < noise] = np.nan

    signal_map2 =  (signal_to_noise(nside, band2))[1]
    m2_map_snr = (signal_to_noise(nside, band2))[0]
    m2 = signal_map2.copy()
    m2[m2_map_snr < noise] = np.nan

    m1 = hp.ma(m1)
    m2 = hp.ma(m2)

    beta_map = (np.log(m2)-np.log(m1))/(np.log(nu[band2])-np.log(nu[band1]))
    
    # Where beta_map is between -2.3 and -1.5 make these areas red, this is freefree are:

    return beta_map, m1, m2

beta_0_1 = maps_masked(0, 1, noise = 20, nside = nside)[0]
beta_2_3 = maps_masked(2, 3, noise = 20, nside = nside)[0]
beta_4_5 = maps_masked(4, 5, noise = 20, nside = nside)[0]
beta_6_7 = maps_masked(6, 7, noise = 20, nside = nside)[0]


lim = 10**-1



def color_params(beta_map):
    #-2.26  manually estimated beta
    freemin = -2.3
    freemax = -2.2
    freefree = np.logical_and(beta_map > freemin, beta_map < freemax)

    beta_map[freefree] = 0

    #-3.34  manually estimated beta
    synchmin = -3.5
    synchmax = -3.1
    synchrotron = np.logical_and(beta_map > synchmin, beta_map < synchmax)

    beta_map[synchrotron] = 150

    #0.808  manually estimated beta
    thermmin = 0.7
    thermmax = 0.9
    thermal_dust = np.logical_and(beta_map > thermmin, beta_map < thermmax)

    beta_map[thermal_dust] = 300

    beta_map[~(freefree | synchrotron | thermal_dust)] = -300

    #hp.gnomview(beta_map, rot = [40], title = f"Nside = {nside}Freefree 0, synchrotron 150, thermal dust 300", unit = "K", format = "%.2g", xsize = 1500, ysize = 500, cmap = 'hot')
    #plt.show()
    #plt.savefig(f"Betamap with colored params for nside = {nside}")

    return beta_map
    

#FIGURING OUT WHERE THE PARAMETERS ARE, AKA WHICH PIXELS


"""Sanity checks"""

position_map = color_params(beta_2_3)

def pixels(color_map, param):
    tol = 0
    #sn_pix = (signal_to_noise(nside, band))[0]
    #np.intersect1d(np.where(color_map > mini)[0], 
    pixels = [np.where(np.isclose(color_map, param, atol=tol))[0]] #cutting 1 and assuming this works
    print(f'pixels {pixels}')
    ooga_boogas = pixels[0] #np.intersect1d(pixels[0], pixels[1]) #np.intersect1d(ooga, boogas)
    return ooga_boogas

idx = pixels(position_map, 0)
print(f'index {idx}')




def spec_plot(index, param):
    m_list = []
    sn_map_list = []

    for i in range(8):
        sn = signal_to_noise(nside, i)[1]
        un = signal_to_noise(nside, i)[2]
        m_list.append(sn[index][0])
        sn_map_list.append(sn)
        
        #plt.scatter(i, sn[index][0])  #picked the first
        plt.errorbar(i, sn[index][0], yerr = un[index][0], fmt = '*')  #picked the first
    plt.title(f'Random pixel containing {param}')
    plt.xlabel('Band')
    plt.ylabel('Signal')
    plt.show()
    return 

run = spec_plot(idx, 'freefree')
# print(f'Shape {np.shape(beta_2_3)}')
# print(f'Shape 2 {np.shape(beta_2_3)[0]}')

# print(f'Betamap {beta_2_3}')

def prob_plot(beta_map):
    hist, bins = np.histogram(beta_map[np.isfinit(beta_map)], bins=np.linspace(-5, 5, 50))
    plt.bar(bins[:-1], hist)
    plt.show()
    return

#test = prob_plot(beta_2_3)
#print(f'Time is {end-start}')

#sn_map = hp.pixelfunc.ma(sn_map_list[0], badval=np.arange(11))



""" make pixel choice based off of masking beta and finding region ang2pix"""


