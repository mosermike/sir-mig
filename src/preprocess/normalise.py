"""
Normalises a data cube
"""

from definitions import *
import numpy as np
from astropy.io import fits
import sys, os
sys.path.append("..")
import sir
import definitions as d
import obs

# TODO SAVE AS PROFILE BINARY FILE
print("TODO AS A BINARY FILE")
sys.exit()

def help():
	print("normalise.py - Normalises a data cube")
	print("Usage: python normalise.py [OPTION]")
	print()
	sir.option("[1. Pos.]","config file")
	sys.exit()

def normalise(conf, stokes= ''):
	"""
	Normalise the data cube
	"""
	if stokes == '':
		stokes = obs.load_data(conf)
		stokes = stokes

	if len(conf['quiet_sun']) > 1:
		print("[STATUS] Normalise data ...")
		# Check if data needs to be normalised
		if np.mean(stokes[:,:,0,:]) < 50:
			print("Is the data already normalised? Abort")
			return

		ll = np.load(os.path.join(conf['path'],conf['waves']))

		if conf['instrument'] in d.ll_lit_norm:
			ll1      = np.argmin(abs(ll-d.ll_lit_norm[conf['instrument']][0]))	 # Find lower limit of wavelength for continuum
			ll2      = np.argmin(abs(ll-d.ll_lit_norm[conf['instrument']][1]))	 # Find upper limit of wavelength for continuum
		else:
			print("[WARN] No instrument defined/Instrument not implemented for Normalisation")
			ll1		= input("Lower Limit of wavelength for continuum in Angstrom: ")
			ll1		= input("Upper Limit of wavelength for continuum in Angstrom: ")

		# Read quiet sun file
		x1	  = conf['quiet_sun'][0]		# Lower limit for region in x
		x2	  = conf['quiet_sun'][1]+1	# Upper limit for region in x
		y1	  = conf['quiet_sun'][2]		# Lower limit for region in y
		y2	  = conf['quiet_sun'][3]+1	# Upper limit for region in y
		
		# Compute continuum intensity in quiet sun region
		Ic  = np.mean(stokes[x1:x2,y1:y2,0,ll1:ll2])  # Average continuum intensity in quiet sun

		# Divide through the mean
		stokes[:,:,0,:] = stokes[:,:,0,:]/Ic
		stokes[:,:,1,:] = stokes[:,:,1,:]/Ic
		stokes[:,:,2,:] = stokes[:,:,2,:]/Ic
		stokes[:,:,3,:] = stokes[:,:,3,:]/Ic

	if conf['fts_file'] == '':
		if ".npy" in conf['cube_inv']:
			np.save(os.path.join(conf['path'],conf['cube_inv']), stokes)
		else:
			# Write the merged data cube
			example = fits.open(os.path.join(conf['path'],conf['cube']))
			header = example[0].header
			hdu = fits.PrimaryHDU(stokes)
			hdu.header = header
			hdu.writeto(os.path.join(conf['path'],conf['cube_inv']),overwrite=True)
	print("-------> Saving data (this might take a while) ...")
	if "npy" in conf['cube']:
		np.save(os.path.join(conf['path'],conf['cube']).replace(".npy",d.end_norm + ".npy"), stokes)
	else:
		# Open the cube
		header = fits.open(os.path.join(conf['path'],conf['cube']))[0].header
		# Write to output
		xxx = fits.PrimaryHDU(stokes)
		xxx.header = header # Copy the header
		xxx.writeto(os.path.join(conf['path'],conf['cube']).replace(".fits",d.end_norm + ".fits"),overwrite=True)

if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])	
	normalise(conf, stokes = '')
	
