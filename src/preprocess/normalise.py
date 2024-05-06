"""
Normalises a data cube
"""

from definitions import *
import numpy as np
import sys
import os
sys.path.append("..")
import sir
import definitions as d
import profile_stk as p

def help():
	print("normalise.py - Normalises a data cube")
	print("Usage: python normalise.py [OPTION]")
	print()
	sir.option("[1. Pos.]","config file")
	sys.exit()

def normalise(conf):
	"""
	Normalise the data cube
	"""
	stokes = p.read_profile(os.path.join(conf["path"], conf["cube"]))

	if len(conf['quiet_sun']) > 1:
		print("[STATUS] Normalise data ...")
		# Check if data needs to be normalised
		if np.mean(stokes.stki) < 10:
			print("Is the data already normalised? Abort")
			return

		ll = stokes.wave

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
		stokes.stki /= Ic
		stokes.stkq /= Ic
		stokes.stku /= Ic
		stokes.stkv /= Ic
	else:
		print("-------> Skipping normalisation")

	print("-------> Saving data (this might take a while) ...")
	if ".bin" in conf['cube']:
		stokes.write(os.path.join(conf['path'],conf['cube']).replace(".bin",d.end_norm))
	else:
		stokes.write(os.path.join(conf['path'],conf['cube']) + d.end_norm)


if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])	
	normalise(conf, stokes = '')
	
