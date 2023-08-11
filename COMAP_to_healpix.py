""" 
                    INSTRUCTIONS

Use batch script process.sh (Mac/Linux ect) or process.bat (Windows) to run this code

- Give your PC permission to run the batch script:
    - Mac/Linux/Unix: Run the line "chmod +x process.sh" in the same folder as the batch script in terminal
    - Windows: Run the line "icacls process.bat /grant Users:F" in the same folder as the batch script in terminal/Power Shell

- Run:    
    - Mac/Linux/Unix: Type "./process.sh" in terminal
    - Windows: Type "process" in terminal/Power Shell

- Make sure python code and batch script is in the same folder as the folder 'Galactic_2021-10-19'

If you don't want to do this, you can replace the filename variable with the file you want to process, and remove './COMAP_to_healpix/' from the savefig and savetxt functions
"""
from astropy.io import fits
import numpy as np
import healpy as hp
import reproject  
import matplotlib.pyplot as plt
import sys

N = 2048 #nside value

"""Collecting data"""
filename = sys.argv[2] 
savename = filename[:-5] #removes the .fits for the filename to be saved in a new format


loadfile = f'./Galactic_2021-10-19/{filename}'


data = fits.open(loadfile) 
hdu = data[3]  #Gives best plots


healpix_hdu = reproject.reproject_to_healpix(hdu, 'Galactic', hdu_in=0, order='nearest-neighbor', nested=False, nside=N)
healpix_map, footprint = healpix_hdu


variance = data[1].data
print(np.shape(variance))
std = np.sqrt(variance) #Variance data from COMAP now represented in standard deviation

"""Overwriting std data onto variance data"""
data[1].data = std
data.writeto('standard_deviation.fits', overwrite=True)
data.close()

"""Write file with standard deviation"""
np.savetxt(f'./COMAP_to_healpix/TXT_files/{savename}_std_dev.txt', std) #Saves the standard deviation data to the folder 'COMAP_to_healpix' created in the batch script used to run this code


"""Create healpix plots in PDF format"""
#SIGNAL MAP:

lim = 10**-1
hp.gnomview(
    healpix_map, 
    title=f"Signal plot for nside = {N}",
    unit="K",
    min=-lim,
    max=lim,
    format="%.2g"
)
hp.graticule()

plt.savefig(f'./COMAP_to_healpix/Plots/{savename}_healpix.pdf', format='pdf', bbox_inches='tight')  #Saves the plot to the folder 'COMAP_to_healpix' created in the batch script used to run this code


#UNCERTAINTY PLOT:


std_fit = fits.open('standard_deviation.fits') 

uncertainty_hdu = reproject.reproject_to_healpix(std_fit, 'Galactic', hdu_in=1, order='nearest-neighbor', nested=False, nside=N)
# print(np.shape(uncertainty_hdu))
uncertainty_map, un_footprint = uncertainty_hdu

# print(np.shape(uncertainty_map))
# print(np.shape(std))

hp.gnomview(
    uncertainty_map,
    title=f"Uncertainty plot for nside = {N}",
    unit="K",
    min= -lim,
    max= lim,
    cmap = 'bwr',
    format="%.2g"
)
hp.graticule()
plt.savefig(f'./COMAP_to_healpix/Plots/{savename}_healpix_standard_deviation_plot.pdf', format='pdf', bbox_inches='tight')  #Saves the plot to the folder 'COMAP_to_healpix' created in the batch script used to run this code


#write maps into fits file 
fits_signal = hp.fitsfunc.write_map(f'./COMAP_to_healpix/Fits_files/Signal_healpix_{savename}.fits', healpix_map, nest=False, dtype=None, fits_IDL=True,  coord='G', partial=False, column_names=None, column_units=None, extra_header=(), overwrite=True)
fits_unc = hp.fitsfunc.write_map(f'./COMAP_to_healpix/Fits_files/Uncertainty_healpix_{savename}.fits', uncertainty_map,  nest=False, dtype=None, fits_IDL=True, coord='G', partial=False, column_names=None, column_units=None, extra_header=(), overwrite=True)

