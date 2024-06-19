"""
Plot several things to get an idea about how the data looks.
"""
import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import sys
sys.path.append(sys.path[0] + "/../src")
import sir
import profile_stk as p
if "-h" in sys.argv:
	print("gris_sketch - Plots various plots for a preliminary analysis for GRIS data")
	print("Usage: python gris_sketch.py [OPTION]")
	print()
	sir.option("1:","data cube name")
	sir.option("2:","path where the figures are saved (optional)")
	print()
	sys.exit()  

# Import library
sir.mpl_library()

filename = sys.argv[1]
savepath = sys.argv[2]

if len(sys.argv) < 2:
     savepath = './'

# Load data
stokes = p.read_profile(filename)
ll1 = np.argmin(abs(stokes.wave-15661.5)) 			 # Find position 15661.5 A

####################################################
# Plot different spectra for different wavelengths #
# to give an overview of the measurement           #
####################################################
wvls = [27,30,40]
imgscale = 7
for l in wvls:
	fig, ax = plt.subplots(figsize=[imgscale * 1.4, imgscale * stokes.ny/stokes.nx])
	ax.imshow(stokes.stki[:,:,l].transpose(), origin = 'lower')
	plt.savefig(savepath + "image"+str(l),bbox_inches='tight', dpi=300)

####################################################
# Plot the spectrum in pixel at a given wavelength #
####################################################
test = 45
fig, ax = plt.subplots(figsize=[imgscale * 1.4, imgscale * stokes.ny/stokes.nx])
plt.imshow(stokes.stki[:,:,test].transpose(), origin='lower')
plt.savefig(savepath + "spectrum.png",bbox_inches='tight', dpi=300)

######################################################
# Plot the mean spectrum in all pixel for comparison #
######################################################
mean = np.mean(stokes.stki[:,:,:],axis=(0,1))

plt.close(fig)

fig, ax = plt.subplots()
plt.plot(stokes.wave,mean/np.amax(mean),label='GRIS')
plt.legend()
plt.tight_layout()
plt.savefig(savepath + "spectrum_comparison.png", bbox_inches='tight', dpi=300)


fig, ax = plt.subplots()
plt.plot(mean/np.amax(mean),label='GRIS')


plt.legend()
plt.savefig(savepath + "spectrum_comparison_p.png", bbox_inches='tight', dpi=300)
plt.show()





