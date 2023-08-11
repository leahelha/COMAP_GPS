from astropy.io import fits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import corner

"""
Run Data_metro.py before running this program.
Data_metro.py gives you the signal and uncerainty data for the specified index you want.
Make sure to specify the correct index in both files

"""

nside = 512

idx = 23646863 # a really good freefree index as a test
A = 0.125 + 0.3
beta = -2.26 + 1 
nu_ = 1
data = np.loadtxt(f'metro_data_index_{idx}.txt')
uncertainty = np.loadtxt(f'metro_uncertainty_index_{idx}.txt')

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

def spec_plot(index, param):
    #Plots spec index plot for one pixel
    band = [26.7, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5]

    data = np.loadtxt(f'metro_data_index_{idx}.txt')
    uncertainty = np.loadtxt(f'metro_uncertainty_index_{idx}.txt')

    for i in range(8):
        plt.errorbar(band[i], data[i], yerr = uncertainty[i], fmt = '*')  #picked the first
    plt.title(f'Random pixel containing {param}')
    plt.xlabel('Band')
    plt.ylabel('Signal')
    #plt.show()
    return 



"""Make metropolis hastings that will calculate"""

#Make TARGET spec index slope
def target(freefree, synchrotron, dust):
    if freefree == 1:
        return -2.26
    
    if synchrotron:
        return -3.34
    
    if dust:
        return 0.808

def log_likelihood(idx, A, beta):
    nu = [26.7, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5]

    nu_ = 30
    log_L = 0
    for i in range(8):
        log_L += (-1/2)*np.sum((data[i]-A*(nu[i]/nu_)**beta)**2/uncertainty[i]**2)

    return log_L



def metro_hastings(idx, A, beta, accepted):

    log_L = log_likelihood(idx, A, beta)

    A2 =  np.random.normal(A, 0.0005, 1)[0] #suggested A with random numbers
    beta2 = np.random.normal(beta, 0.01, 1)[0]
    #print(f'B2 is {beta2}')
    log_L2 = log_likelihood(idx, A2, beta2)

    epsilon = np.random.uniform(0,1)  
    acceptance = np.min((1, np.exp(log_L2-log_L)))

   
    if acceptance > epsilon:
        amplitude = A2
        spectral_index = beta2
        accepted += 1
        
    else:
        amplitude = A
        spectral_index = beta

    A_final = amplitude
    beta_final = spectral_index


    return A_final, accepted, beta_final



def cycles(steps, idx, A, beta, accepted):
    A_list = []
    beta_list = []
    

    for i in range(steps):
        A, accepted, B = metro_hastings(idx, A, beta, accepted)

        A_list.append(A)
        beta_list.append(B)
        beta = B

    acceptance_rate = accepted/steps
   


    return A_list, acceptance_rate, beta_list




plt.figure()

run = 8001

beta_list = np.linspace(-1, -3, 30)


spec_plot = spec_plot(idx, 'freefree')
metropolis_A = cycles(run, idx, A, beta, 0)[0]
metropolis_beta = cycles(run, idx, A, beta, 0)[2]
print(f'BETA METROPOLIS {metropolis_beta[7990:]}')


acceptance_rate = cycles(run, idx, A, beta, 0)[1]
print(f'Acceptance rate is {acceptance_rate*100}%')
#print(metropolis_A)


band = [26.7, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5]
cycle100 = [metropolis_A[-3000]*(band[i]/nu_)**metropolis_beta[-3000] for i in range(8)]
cycle200 = [metropolis_A[-1000]*(band[i]/nu_)**metropolis_beta[-1000] for i in range(8)]
cycle300 = [metropolis_A[-1]*(band[i]/nu_)**metropolis_beta[-1] for i in range(8)]

plt.plot(band, cycle100, label = '1')
plt.plot(band, cycle200, label = '2')
plt.plot(band, cycle300, label = '3')

plt.legend()

#print(metropolis_A[-1])

plt.figure()
plt.plot(metropolis_A)

plt.figure()
plt.plot(metropolis_beta)


plt.figure()
plt.hist(metropolis_A[2500:], bins = 50)


plt.figure()
plt.hist(metropolis_beta[2500:], bins = 50)



fig = corner.corner(np.array([metropolis_A[2500:], metropolis_beta[2500:]]).T)
plt.show()




#Shows rough estimate areas of synchrotron, freefree or thermal dust by comparing two bands
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

# beta_0_1 = maps_masked(0, 1, noise = 20, nside = nside)[0]
# beta_2_3 = maps_masked(2, 3, noise = 20, nside = nside)[0]
# beta_4_5 = maps_masked(4, 5, noise = 20, nside = nside)[0]
# beta_6_7 = maps_masked(6, 7, noise = 20, nside = nside)[0]
