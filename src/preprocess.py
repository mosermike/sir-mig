
"""

Preprocessing
=============

These module provides functions to preprocess data such as merging, normalisation and correction for a spectral veil for groundbased telescopes.

Merging Data
------------

Merging data downloaded from SDC to a data cube or Hinode data. It puts the 3D arrays into a 4D array to be easier
accessible. In total, four different measurements can be put into 4 data cubes at the same time. For more measurements
in the same folder, the script must be adapted.

Normalise Data
--------------

Normalises a data cube by computing the average intensity from a quiet sun region

Spectral veil correction
------------------------

Convolve the FTS data with a Gaussian and a constant value for the spectral veil using
	$$I = (1 - \\nu) [I_{FTS} * g(\\lambda, \\sigma)] + \\nu I_c.$$

It computes the best fit parameter sigma and nu. Afterwards, it corrects the fits data.
It writes the results into single files and into a corrected data cube.

"""

import numpy as np
from astropy.io import fits
import os
import sys
import matplotlib.pyplot as plt
from os.path import exists
from astropy.convolution import convolve_fft
from tqdm import tqdm

import sir
import definitions as d
import profile_stk as p

def merge(dir : str, ending : str, instrument : str, path = "./", shift = "0", save = False):
	"""
	
	Merges data to a cube

	Parameters
	----------
	dir : str
		Directory of the fits
	ending : str
		Ending of the dataset (for GRIS data)
	instrument : str
		Used instrument (Implemented Instruments: GRIS or Hinode)
	path : str
		Path where files are written
	shift : str
		Shift the wavelengths by this value in mA
	save : bool
		Save the merged data cube
	Returns
	-------
	merge : profile_stk
		Class with the merged profiles

	"""
	if instrument not in d.instruments:
		print(f"-------> Instrument {instrument} is not known. Consider adding it to the definitions.py script.")
		print(f"-------> The code assumes that all files in {dir} ending with '.fits' belong to the data.")
		print(f"-------> The code will try to make it work and will ask you for some inputs.")
		print(f"-------> Consider adapting the script yourself or make an issue in the corresponding git website.")

	output = os.path.join(path, d.cube_merge)  # where it is saved

	# Files in the folders from different measurements
	filenames = []

	if instrument == 'GRIS':
		print(f"-------> Use files containing _{ending}_0 in the file name")

	# Add all fits files in a list
	for File in sorted(os.listdir(dir)):
		if File.endswith(".fits"):
			if ("_" + ending + "_0" in File and instrument == 'GRIS') or (instrument != 'GRIS'):
				filenames.append(os.path.join(dir,File))

	####################
	# Create data cube #
	####################

	example = fits.open(filenames[0])
	stokes = example[0].data
	header = example[0].header

	#########################################
	# 	Load data and store it in an array	#
	#########################################
	# Get no. of wavelengths and pixels ; Reshape the array in a different order
	# ADAPT for another instrument here for the structure
	if instrument == 'Hinode':
		nw = stokes.shape[2]
		nx = len(filenames)
		ny = stokes.shape[1]
	elif instrument == 'GRIS':
		nw = stokes.shape[1]
		nx = len(filenames)
		ny = stokes.shape[2]
	else:
		temp = input("Specify which entry corresponds to (nw and ny) as a list: [for GRIS: 1,2] ").split("'")
		print("-------> The code assumes that the number of files corresponds to the x direction.")
		nw = stokes.shape[int(temp[0])]
		nx = len(filenames)
		ny = stokes.shape[int(temp[1])]

	print('# Pixels in x  = ', nx)
	print('# Pixels in y  = ', ny)
	print('# Wavelengths  = ', nw)

	data = np.empty((nx, ny, 4, nw))  # shape 4, wavelength, y, 1

	# Create data cube
	for i in tqdm(range(len(filenames))):
		# Load data from one file
		example = fits.open(filenames[i])

		stokes = example[0].data
		if instrument == d.instruments[1]:
			spbshft = header['SPBSHFT']
			stokes = example[0].data.astype(np.float32)
			stokes = np.transpose(stokes, axes=(1, 0, 2))

			# Correction
			if spbshft > 0:
				stokes[:, 0, :] *= 2
				negative = np.where(stokes[:, 0, :] < 0.0)
				stokes[:, 0, :][negative] += 65536.

			if spbshft > 1:
				stokes[:, 3, :] *= 2
			if spbshft > 2:
				stokes[:, 1, :] *= 2
				stokes[:, 2, :] *= 2

			data[i] = stokes[:, :, :]

		elif instrument == d.instruments[0]:
			stokes = np.transpose(stokes, axes=(3, 2, 0, 1))  # Structure of GRIS stokes, nw, ny, nx
			data[i] = stokes[0, :, :, :]
		else:  # ADAPT for another instrument here for the structure
			print("-------> It assumes the structure as GRIS if it has dim 4.")
			print("-------> It assumes the structure as Hinode if it has dim. 3.")
			if len(stokes.shape) == 4:
				stokes = np.transpose(stokes, axes=(3, 2, 0, 1))  # wanted structure nx, ny, stokes, nw
			elif len(stokes.shape) == 3:
				stokes = np.transpose(stokes, axes=(1, 0, 2))  # wanted structure nx, ny, stokes, nw
			else:
				print("        Manual adaptation necessary!")
				return

			data[i] = stokes[0, :, :, :]

	data = data.reshape(nx, ny, 4, nw)  # shape x, y, values, wavelengths
	print()
	##############################
	# Determine real wavelengths #
	##############################
	if instrument == 'GRIS':
		print("[STATUS] Create wavelength arrays assuming the header as for GRIS data")
		ll_a = header["CRVAL3"]
		ll_b = header["CDELT3"]
		llambda = ll_a + ll_b * np.arange(0, header["NAXIS3"])  # Measured wavelength for each pixel
	elif instrument == 'Hinode':
		#llambda = np.linspace(6300.87730065, 6303.25996109, 112)
		llambda = np.float32((np.arange(112)-np.float(header['CRPIX1']))*abs(header["CDELT1"]) + header["CRVAL1"])
	else:
		print(f"Instrument {instrument} not known. Define wavelength grid manually:")
		mins = input("Lower wavelength in A: ")
		maxs = input("Upper wavelength in A: ")
		values = input("Number of spectral points: ")
		llambda = np.linspace(float(mins), float(maxs), int(values))

	if shift!= '0' or shift != "":
		print(f"Wavelength Grid shifted by {shift} mA.")
		llambda = llambda + float(shift) * 1e-3
	llambda = np.float32(llambda)
	print(f"Wavelength step is {'%.3f' % ((llambda[1]-llambda[0])*1e3)} mA.")

	print("-------> Assign data ...")
	pro = p.profile_stk(nx=data.shape[0],ny=data.shape[1],nw=data.shape[3])
	pro.wave = llambda
	pro.stki = data[:,:,0,:]
	pro.stkq = data[:,:,1,:]
	pro.stku = data[:,:,2,:]
	pro.stkv = data[:,:,3,:]
	pro.load = True # Data is loaded (no warning when saved)

	if save:
		print("-------> Saving data (this might take a while) ...")
		pro.write(output)
		print("Saved as \"%s\"" % output)

	# Save header
	filename = os.path.join(path, d.header_infos)
	with open(filename, 'w') as f:
		f.write(f"instrument={instrument}\n")
		for card in header.cards:
			if card.keyword != 'COMMENT':
				f.write(f"{card.keyword}={card.value}\n")
		f.write(f"SHIFT={shift}\n")
		if instrument == 'GRIS':
			f.write(f"END={ending}") # Ending of the used files
		#elif conf['instrument'] == 'Hinode':

	return pro

def merge_conf(dir : str, conf : dict):
	"""
	Merges data to a cube with the config file

	Parameters
	----------
	dir : str
		Directory of the fits
	config : dict
		Dictionary with all the information from the config file
	
	Return
	------
	merge : profile_stk
		Class with the merged profiles
	
		
	"""
	return merge(dir, conf["ending"], conf["instrument"], conf["path"], conf["shift_wave"], conf["save_cube"] == "1")

#################
#	NORMALISE	#
#################


def normalise(pro : p.profile_stk, instrument : str, quiet_sun : list, path = "./", save=False) -> p.profile_stk:
	"""
	Normalise the data cube by the given quiet sun range

	Parameters
	----------
	pro : profile_stk
		Data cube with the Stokes Profiles
	instrument : str
		Used instrument (Implemented Instruments: GRIS or Hinode)
	quiet_sun : list
		List with the pixel for the quiet sun as xmin,xmax,ymin,ymax
	path : str
		Path where files are written
	save : bool
		Save the merged data cube

	Return
	------
	normalise : profile_stk
		Normalised data cube
	"""
	if len(quiet_sun) > 1:
		print("[STATUS] Normalise data ...")
		# Check if data needs to be normalised
		if np.mean(pro.stki) < 10:
			print("Is the data already normalised? Abort")
			return

		ll = np.copy(pro.wave)

		if instrument in d.ll_lit_norm:
			ll1      = np.argmin(abs(ll-d.ll_lit_norm[instrument][0]))	 # Find lower limit of wavelength for continuum
			ll2      = np.argmin(abs(ll-d.ll_lit_norm[instrument][1]))	 # Find upper limit of wavelength for continuum
		else:
			print("[WARN] No instrument defined/Instrument not implemented for Normalisation")
			ll1		= input("Lower Limit of wavelength for continuum in Angstrom: ")
			ll1		= input("Upper Limit of wavelength for continuum in Angstrom: ")

		# Read quiet sun file
		x1	  = quiet_sun[0]	# Lower limit for region in x
		x2	  = quiet_sun[1]+1	# Upper limit for region in x
		y1	  = quiet_sun[2]	# Lower limit for region in y
		y2	  = quiet_sun[3]+1	# Upper limit for region in y

		# Compute continuum intensity in quiet sun region
		Ic  = np.mean(pro.stki[x1:x2,y1:y2,ll1:ll2])  # Average continuum intensity in quiet sun

		# Divide through the mean
		pro1 = pro.copy()
		pro1.stki /= Ic
		pro1.stkq /= Ic
		pro1.stku /= Ic
		pro1.stkv /= Ic
	else:
		print("-------> Skipping normalisation")

	if save:
		print("-------> Saving data (this might take a while) ...")
		pro1.write(os.path.join(path,d.cube_norm))

	return pro1

def normalise_conf(pro : p.profile_stk, conf : dict) -> p.profile_stk:
	"""
	Normalise the data cube by the given quiet sun range with the config file

	Parameters
	----------
	pro : profile_stk
		Data cube with the Stokes Profiles
	config : dict
		Dictionary with all the information from the config file

	Return
	------
	normalise : profile_stk
		Normalised data cube
	

	"""
	return normalise(pro, conf["instrument"], conf["quiet_sun"], conf["path"], conf["save_cube"] == "1")

#################################
#	SPECTRAL VEIL CORRECTION	#
#################################

def argmin(x : np.array):
	r"""
	Find the argument of the minimum in an multi-dimensional
	array. It seems to be faster than any numpy-built-in 
	function as shown in https://stackoverflow.com/questions/30180241/numpy-get-the-column-and-row-index-of-the-minimum-value-of-a-2d-array.

	Parameters
	----------
	x : $m$x$n$ numpy array
		Numpy array with the dataq
	
	Returns
	-------
	argmin : int
		First index in dim. $m$
	argmin : int
		Second index in dim. $n$

	"""
	k = x.argmin()
	ncol = x.shape[1]
	return int(k/ncol), k%ncol

#################################################################################################3

def chi2(y_fit : np.array, y_obs : np.array):
     r"""
     Computes the merit-function $\chi^2$. The used equation is $$\chi^2 = \sum_i y_{\text{fit}} - y_{\text{obs}}.$$
     
     Parameters
     ----------
     y_obs : numpy array
          Array with the observed values
     y_fit : numpy array
          Array with computed values, e.g. from a fit

     Returns
     -------
     chi2 : float
          $\chi^2$-value
     
     """
     return np.sum((y_fit-y_obs)**2)

#################################################################################################3

def gaussian(x : np.array, mean = 0, sigma = 1, norm = True):
	r"""
	Computes the value of a Gaussian at position x.
	
	Parameters
	----------
	x : float
		Wavelength
	mean : float
		Mean of the Gaussian
	sigma : float
		Width of the Gaussian
	norm : bool, optional
		Normalise by $\frac{1}{\sqrt(2\pi\sigma^2]}$. Default is True

	Returns
	-------
	gaussian : float
		Value in the Gaussian function at position x
	
	"""
	g = np.exp(-(x - mean)**2 / (2 * sigma**2))
	if norm:
		g = g / (sigma * np.sqrt(2 * np.pi)) # Normalised so that int g = 1
	return g

#################################################################################################3

def convective_blueshift(ll : np.array,I : np.ndarray):
	"""
	Determines the shifted $\\lambda_0$ due to convective blueshift.

	Parameters
	----------
	ll : np.array
		Wavelengths in a given range
	I  : numpy.ndarray
		Stokes I in the given range

	Returns
	-------
	convective_blueshift : int	
		Index of the mean shifted spectral core
	convective_blueshift : float
		Wavelength of the mean shifted spectral core

	"""
	Min = []
	for i in range(I.shape[0]):
		for j in range(I.shape[1]):
			Min.append(np.argmin(I[i,j,:]))

	# Interpolate to find minimum
	l0 = np.interp(np.mean(Min), np.arange(0,len(ll),1), ll)

	return int(np.mean(Min)), l0

#################################################################################################3

def optimise_chi(nu : np.array, sigma : np.array, I : np.ndarray, I_obs : np.ndarray):
	"""
	Optimises two parameters by computing the chi2 merit function and finding
	the minima of the two parameters. The steps are the following
	
	 1. Compute the chi2 values
	 2. Compute the discrete minima in each parameter
	 3. Perform two polynomial fits to find the optimal minimum
	 4. Compute the uncertainty of the fit
	
	Parameters
	----------
	nu : numpy array
		Array containing the nus which were used
	sigma :  numpy array
		Array containing the sigmas which were used
	I :  numpy ndarray
		Multidimensional array with the following format: I[nu, sigma, wavelength]. This are the fitted/computed data to be compared to the observed one
	I_obs :  numpy array
		Array containing the observed intensity
	
	Returns
	-------
	optimise_chi : float
		Optimised sigma
	optimise_chi : float
		Uncertainty of the optimised sigma
	optimise_chi : float
		Optimised nu
	optimise_chi : float
		Uncertainty of the optimised nu
	optimise_chi : ndarray
		Array containing all the chi2 for all possible combinations of nu and sigma
	
	"""
	chis = np.zeros(shape=(len(nu), len(sigma)))
	for n in range(len(nu)):
		for s in range(len(sigma)):
			chis[n,s] = np.sum((I[n,s,:]-I_obs)**2)

	# Find the simple minima
	chi_min = argmin(chis)
	nu_min = nu[chi_min[0]]
	sigma_min = sigma[chi_min[1]]

	if chi_min[0]+5 > len(nu) or chi_min[0] < 4:
		raise Exception("[optimise_chi] The found minimum (" + str(chi_min[0]) + "," + str(nu_min) + ") is at the border of the selected ranges. Did you select a QS region?")
	if chi_min[1]+5 > len(sigma) or chi_min[1] < 4:
		raise Exception("[optimise_chi] The found minimum (" + str(chi_min[1]) + "," + str(sigma_min) + ") is at the border of the selected ranges. Did you select a QS region?")

	
	##################
	# Optimise sigma #
	##################
	pol = np.polyfit(sigma[chi_min[1]-4 : chi_min[1]+4],
				  chis[chi_min[0], chi_min[1]-4 : chi_min[1]+4], 2, cov='unscale')
	pcov = pol[1]
	pol = pol[0]
	# Now, the pixel position where the parabola is minimum is -b/(2*c)
	sigma_min = -pol[1] / (2 * pol[0])
	usigma_min = np.sqrt(
					pcov[1][1] * (-1/ (2 * pol[0]))**2 + 
					pcov[0][0] * (pol[1] / (2 * (pol[0])**2))**2
					)

	###############
	# Optimise nu #
	###############
	pol = np.polyfit(nu[chi_min[0]-4 : chi_min[0]+4], chis[chi_min[0]-4 : chi_min[0]+4, chi_min[1]], 2, cov='unscale')
	pcov = pol[1] # Covariance matrix
	pol = pol[0]  # Polynomials fit parameters
	# Now, the pixel position where the parabola is minimum is -b/(2*c)
	nu_min = -pol[1] / (2 * pol[0])
	unu_min = np.sqrt(
					pcov[1][1] * (-1/ (2 * pol[0]))**2 + 
					pcov[0][0] * (pol[1] / (2 * (pol[0])**2))**2
					)
	return sigma_min, usigma_min, nu_min, unu_min, chis

#################################################################################################3

def vac_to_air(wavelength : float, method = "Ciddor1996") -> float:
	r"""
	Computes the wavelength from vacuum to air by using (Ciddor, 1996) eq. 1.

	Note that this equation is only valid for wavelengths between $2300 A$ and $16900 A$ and
	it should be better in the infrared than Edlen.

	**Note that the CO2 concentration is here assumed to be 450 ppm but it is most likely not constant in general.**
	
	Parameters
	----------
	wavelength : float
		Wavelength in Angstrom to be converted
	method : str, optional
		Determines which method is used. Only the method "Ciddor1996" is implemented so far. Default: "Ciddor1996"


	Returns
	-------
	vac_to_air : float
		Corresponding wavelength in air in Angstrom

	"""
	if method == "Ciddor1996":
		sigma2 = (1e4 / wavelength)**2 # squared wavenumber in mum^-1
		refr    = (1 + 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2))
		wave_air = wavelength/refr
	else:
		print(f"[vac_to_air] Method {method} not defined.")
	return wave_air


#################################################################################################3
def correct_spectral_veil_conf(pro : p.profile_stk, conf : dict) -> p.profile_stk:
	""" 
	Correct the spectral veil in the data with a config file. This function calls the following functions:
	 - argmin()
	 - chi2()
	 - gaussian()
	 - convective_blueshift()
	 - optimise_chi()
	 - vac_to_air()

	This function convolves the FTS data with a Gaussian and a constant value for the spectral veil using
	$$I = (1 - \\nu) [I_{FTS} * g(\\lambda, \\sigma)] + \\nu I_c.$$

	Parameters
	----------
	pro : profile_stk
		Class with the normalised Stokes Profiles
	config : dict
		Dictionary with all the information from the config file

	Return
	------
	correct_spectral_veil_conf : profile_stk
		Corrected profiles

	"""
	return correct_spectral_veil(pro, conf["instrument"], conf["fts_file"], conf["quiet_sun"], conf["cube"], conf["path"])
	
def correct_spectral_veil(pro : p.profile_stk, instrument : str, fts_file : str, quiet_sun : list, cube : str, path : str):
	"""
	Correct the spectral veil in the data. This function calls the following functions:
	 - argmin()
	 - chi2()
	 - gaussian()
	 - convective_blueshift()
	 - optimise_chi()
	 - vac_to_air()

	This function convolves the FTS data with a Gaussian and a constant value for the spectral veil using
	$$I = (1 - \\nu) [I_{FTS} * g(\\lambda, \\sigma)] + \\nu I_c.$$

	Parameters
	----------
	config : dict
		Dictionary with all the information from the config file
	pro : profile_stk
		Class with the normalised Stokes Profiles
	instrument : str
		Used instrument (Implemented Instruments: GRIS or Hinode)
	fts_file : str
		File to the FTS file
	quiet_sun : list
		List with the pixel for the quiet sun as xmin,xmax,ymin,ymax
	cube : str
		Name of the stored data cube
	path : str
		Path where the files are stored

	Return
	------
	correct_spectral_veil : profile_stk
		Profiles with the corrected spectal veil

	"""
	
	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = os.path.join(path,sys.argv[sys.argv.index("-save")+1])
		if not exists(savepath[:savepath.rfind('/')]):
			os.mkdir(savepath[:savepath.rfind('/')])

	# Do not do anything
	if fts_file == '':
		print("-------> No spectral veil correction, write the merged or normalised data cube")
		#if not exists(os.path.join(conf["path"], conf["cube_inv"])):
			#if len(conf['quiet_sun']) < 2 and exists(os.path.join(conf["path"], conf["cube"])):
			#	print("[STATUS] Copying the non-preprocessed data cube to the file used for the inversion.")
			#	shutil.copy(os.path.join(conf['path'],conf['cube']), os.path.join(conf['path'],conf['cube_inv'])) 
			#	print("[Preprocess] Preprocess data is done. Consider changing from 'preprocess : 1' to 'preprocess : 0'!")
			#elif exists(os.path.join(conf["path"], conf["cube"]) + d.end_norm):
			#	print("[STATUS] Copying the created normalised data cube to the file used for the inversion.")
			#	shutil.copy(os.path.join(conf['path'],conf['cube']).replace(".bin","") + d.end_norm, os.path.join(conf['path'],conf['cube_inv'])) 
			#	print("[Preprocess] Preprocess data is done. Consider changing from 'preprocess : 1' to 'preprocess : 0'!")
			#else:
			#	print(f"[ERROR] File {conf['cube_inv']} is not created. Are there data in the selected path?")
		pro.write(os.path.join(path, cube))
		
		del pro
		return
	
	print("[STATUS] Correct spectral veil ...")
	if exists(os.path.join(path,cube)):
		temp = ''
		print("[WARN] The data cube used for the inversion already exists and spectral veil correction is selected.")
		while temp != 'y' and temp != 'n':
			temp = input("Do you want to overwrite it and continue? [y/n] ")
		if temp == 'n':
			print("Abort (Consider changing the config file)")
			return

	if(instrument != d.instruments[0]):
		print("| Note that the script was only tested for GRIS data in the infrared range.")
		print("| Take a look at the plotted plots and if they do not look good, consider changing the ranges or the borders")
	
	if abs(pro.wave[0] - 15600) > 100 and instrument == d.instruments[0]: # GRIS instrument
		print("| Note that the script was only tested for GRIS data in the infrared range (15600 A).")
		print("| Take a look at the plotted plots and if they do not look good, consider changing the ranges or the borders")
		print("| Eventually add another instrument named e.g. 'GRIS1' for the different wavelength range")

	##########################
	#	Plot settings		#
	########################## 
	sir.mpl_library()

	#################################################################
	#						    DEFINE VARIABLES					#
	#################################################################
	b = 1.5 # For the used range in the plotted spectrum in A
	sigma = np.arange(d.sigma_range[0],d.sigma_range[1]) * 1e-3 # from mA to A
	nu    = np.arange(d.nu_range[0],d.nu_range[1],1) * 1e-2 # from percent to 1

	#################################################################
	#						    INITIALIZATION						#
	#################################################################
	stokes = pro
			
	filename_fts  = fts_file

	# Read quiet sun file
	x1	  = quiet_sun[0]		# Lower limit for region in x
	x2	  = quiet_sun[1]+1	# Upper limit for region in x
	y1	  = quiet_sun[2]		# Lower limit for region in y
	y2	  = quiet_sun[3]+1	# Upper limit for region in y

	# FTS
	data = np.loadtxt(filename_fts)
	wn = np.flip(data[:,0]) # Wavenumber in cm^-1
	i1 = np.flip(data[:,1]) # Intensity corrected


	# GRIS
	ll_gris = np.copy(stokes.wave)

	####################################################################
	#						   SHIFT OF SPECTRUM					   #
	####################################################################
	print('-------> Shift spectrum ...')
	
	#########
	#  FTS  #
	#########
	# Literature wavelength in air where the line core is expected in air
	if instrument in d.ll_lit:
		ll_lit = d.ll_lit[instrument]
	else:
		ll_lit = float(input("Used instrument not defined in definitions. Type the literature wavelength in air, used for the spectral veil correction: "))
	#ll_lit = 15648.515 # From A new Multiplet table for Fe I; imported from definitions

	# Compute wavelength in vacuum and in air
	ll = 1e8 / wn # in wavelength in A
	ll_ciddor = vac_to_air(ll, "Ciddor1996") # Ciddor, 1996 method

	#################
	#    Data		 #
	#################
	# Compute average intensity in quiet sun in data and shift it to lambda_0
	i_gris = np.mean(stokes.stki[x1:x2,y1:y2,:], axis=(0,1))

	# Correct convective blueshift / Shift in GRIS data
	ll1_lit = np.argmin(abs(ll_gris-ll_lit))		 # Find position at ll_lit defined in definitions.py
	border = 40								         # Find minima around this value of the literature value, if this peak is not prominent, than put a smaller number
	ll1, l0 = convective_blueshift(ll_gris[ll1_lit-border:ll1_lit+border], stokes.stki[x1:x2,y1:y2,ll1_lit-border:ll1_lit+border]) # Determine the minima position and value of the shifted GRIS spectrum
	ll_gris += ll_lit-l0 # Shift the wavelength

	######################################################################################
	#					Changing to RELATIVE WAVELENGTHS					  #
	######################################################################################
	# Change lambda to be relative to the literature value
	ll	     -= ll_lit
	ll_ciddor -= ll_lit
	ll_gris   -= ll_lit

	######################################################################################
	#					 MOVING ALL TO ZERO RIGOROUSLY					    #
	######################################################################################
	# Moving intensities to zero so that they really overlap due to any other reasons
	j = 20
	# ciddor intensities
	ll_min = np.argmin(abs(ll_ciddor))
	i_range = i1[ll_min - j : ll_min + j] # Intensity range to check for minima
	ll_ciddor -= ll_ciddor[np.argmin(i_range) + ll_min - j]
	# gris intensities
	ll_min = np.argmin(abs(ll_gris))
	i_range = i_gris[ll_min - j : ll_min + j] # Intensity range to check for minima
	ll_gris -= ll_gris[np.argmin(i_range) + ll_min - j]

	######################################################################################
	#							  CONVOLUTION							 #
	######################################################################################
	# By using Ciddor for the shift
	print("[STATUS] Start convolution ...")

	# Shorten data to the peak from -(-1.5) to 1.5
	ll_conv	= ll_ciddor[np.argmin(abs(ll_ciddor + b)) : np.argmin(abs(ll_ciddor - b))]
	i_conv	= i1	   [np.argmin(abs(ll_ciddor + b)) : np.argmin(abs(ll_ciddor - b))]

	# Interpolate the GRIS data to be comparable with the FTS data
	I_obs = np.interp(ll_conv, ll_gris, i_gris)

	# Shift to zero as it seems to be still shifted
	ll_conv -= ll_conv[np.argmin(i_conv)]
	ll_conv_g = ll_conv - ll_conv[np.argmin(I_obs)]

	# Create Gaussians considering the veil and the width
	Ic = i_conv[0] # Depends on the range of the shortened intensities from the FTS
	Ic = 1. # Seems to be better fitting the wings. In theory, I_c should be 1.
	gs = np.zeros(shape=(len(nu), len(sigma), len(ll_conv)))
	for s in range(len(sigma)):
			gs[:,s,:] = gaussian(ll_conv, sigma = sigma[s])

	# Convolve data by using the fast fourier algorithm
	I_conv = np.zeros(shape=(len(nu), len(sigma), len(ll_conv)))
	for n in range(len(nu)):
		for s in range(len(sigma)):
			I_conv[n,s,:] = (1-nu[n])*convolve_fft(i_conv, gs[n,s], fill_value = Ic) + nu[n]*Ic

	######################################################################################
	#						    CHI2 MERIT FUNCTION						 #
	######################################################################################
	# Compute the chi2-merit function to compute the best fit
	sigma_min, usigma_min, nu_min, unu_min, chis = optimise_chi(nu, sigma, I_conv, I_obs)

	print("υ_min = (%.2f ± %.2f) %%" % (nu_min/1e-2, unu_min / 1e-2))
	print("σ_min = (%.2f ± %.2f) mÅ" % (sigma_min/1e-3, usigma_min / 1e-3))

	# Save the parameters for later use
	np.save(os.path.join(path,d.veil_parameters), [[nu_min, unu_min],[sigma_min,usigma_min]]) # sigma in Angstrom

	# Compute the convolved values with the best fit
	I_conv_best = (1-nu_min)*convolve_fft(i_conv, gaussian(ll_conv, sigma = sigma_min), fill_value = Ic) + nu_min*Ic

	# Shift convolved data so that they are also centered
	LL_conv = ll_conv - ll_conv[np.argmin(I_conv_best)]

	
	
	chi_min = chi2(I_conv_best,I_obs)

	######################################################################################
	#						  Plot the PARAMETER SPACE					   #
	######################################################################################
	fig, ax = plt.subplots()
	vmin = round(np.log10(chi_min)-0.05,1)
	vmax = round(np.log10(chi_min)-0.05,1) + 0.8
	c = ax.contourf(sigma*1e3, nu*100, np.log10(chis), levels = np.linspace(vmin, vmax, 40), cmap = 'ocean', vmin = vmin, vmax = vmax)
	ax.set_ylabel(r"$\nu$ $[\%]$")
	ax.set_xlabel(r"$\sigma$ [m\AA]")
	cbar = fig.colorbar(c,ticks=np.arange(vmin, vmax+0.2, 0.2), location = 'top')
	cbar.set_label(label=r"$\log(\chi^2)$", loc = 'center', labelpad = 15)
	plt.text((d.sigma_range[1]-d.sigma_range[0])*120/140, (d.nu_range[1]-d.nu_range[0])*25/31,
		  		r"$\hat{\sigma} =$ " + "(%.1f ± %.1f)"
					% (round(sigma_min/1e-3,1), round(usigma_min / 1e-3 + 0.05,1)) + r" m\AA"
				+ '\n' +
				r"$\hat{\nu} =$ " + "(%.1f ± %.1f) "
					% (round(nu_min/1e-2,1), round(unu_min / 1e-2 + 0.05,1))  + r"$\%$"
				+ '\n' +
				r"$\hat{\chi}^2 =$ " + "%.3f"
					% (chi_min), fontsize=12)
	plt.savefig(os.path.join(path,savepath + "veil_parameter_space"))


	######################################################################################
	#					    Plot the best convolved intensity				  #
	######################################################################################
	fig, ax = plt.subplots()
	ax.plot(ll_conv, i_conv, "x", label = r'$I^{\mathrm{FTS}}$')
	if instrument == 'GRIS':
		ax.plot(ll_conv_g, I_obs, '+', label = r'$I_{\mathrm{QS}}^{\mathrm{GRIS}}$')
	else:
		ax.plot(ll_conv_g, I_obs, '+', label = r'$I_{\mathrm{QS}}^{\mathrm{data}}$')
	ax.plot(LL_conv, I_conv_best, label = r"$I^{\mathrm{fit}}$")

	ax.set_xlabel(r"$\Delta \lambda^{\mathrm{air}}$" + f" - {d.ll_lit[instrument]}" + r" \AA")
	ax.set_ylabel(r"Intensity $I/I_c$")
	ax.set_xlim(-b,b)

	ax.legend(fontsize=15)

	plt.savefig(os.path.join(path, savepath + "veil_intensity"))

	######################################################################################
	#					START OF THE ACTUAL CORRECTION					    #
	######################################################################################
	print('[STATUS] Correct data ...')
	pro1 = stokes.copy()
	# Correct for spectral veil
	pro1.veil_correction(nu_min)
	

	print("-------> Saving data (this might take a while) ...")
	pro1.write(os.path.join(path,cube))

	print("[Preprocess] Preprocess data is done. Consider changing from 'preprocess : 1' to 'preprocess : 0'!")
	
	return pro1
