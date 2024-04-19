"""
Plot several things to get an idea about how the data looks.
"""
import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import sys
import os
sys.path.append(sys.path[0] + "/..")
import os, sys, sir, obs
import definitions as d
from astropy.io import fits
	
# Import library
dirname = os.path.split(os.path.abspath(__file__))[0]
if exists(dirname + '/../mml.mplstyle'):
	plt.style.use(dirname + '/../mml.mplstyle')
elif "mml" in plt.style.available:
	plt.style.use('mml')




filename = sys.argv[1]
savepath = sys.argv[2]
waves = np.load(sys.argv[3])

if len(sys.argv) < 2:
     savepath = './'

# Load data
stokes = obs.load_data({}, filename = filename, add_path = False)
ll1 = np.argmin(abs(waves-15661.5)) 			 # Find position 15661.5 A

####################################################
# Plot different spectra for different wavelengths #
# to give an overview of the measurement           #
####################################################
wvls = [27,30,40]
imgscale = 7
for l in wvls:
	fig, ax = plt.subplots(figsize=[imgscale * 1.4, imgscale * stokes.shape[1]/stokes.shape[0]])
	ax.imshow(stokes[:,:,0,l].transpose(), origin = 'lower')
	plt.savefig(savepath + "image"+str(l),bbox_inches='tight', dpi=300)

####################################################
# Plot the spectrum in pixel at a given wavelength #
####################################################
test = 45
fig, ax = plt.subplots(figsize=[imgscale * 1.4, imgscale * stokes.shape[1]/stokes.shape[0]])
plt.imshow(stokes[:,:,0,test].transpose(), origin='lower')
plt.savefig(savepath + "spectrum.png",bbox_inches='tight', dpi=300)

######################################################
# Plot the mean spectrum in all pixel for comparison #
######################################################
mean = np.mean(stokes[:,:,0,:],axis=(0,1))

plt.close(fig)

fig, ax = plt.subplots()
plt.plot(waves,mean/np.amax(mean),label='GRIS')
plt.legend()
plt.tight_layout()
plt.savefig(savepath + "spectrum_comparison.png", bbox_inches='tight', dpi=300)


fig, ax = plt.subplots()
plt.plot(mean/np.amax(mean),label='GRIS')
#plt.plot([453,453],[0.6,1.0])
#plt.plot([501,502],[0.6,1.0], '--')
#plt.plot([550,550],[0.6,1.0], '-.')

plt.legend()
plt.savefig(savepath + "spectrum_comparison_p.png", bbox_inches='tight', dpi=300)
plt.show()





