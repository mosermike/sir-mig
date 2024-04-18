"""
Convolve the FTS data with a Gaussian and a constant value for the spectral veil using
	I = (1 - nu) [I_FTS * g(lambda, sigma)] + nu Ic.
It computes the best fit parameter sigma and nu. Afterwards, it corrects the fits data.
It writes the results into single files and into a corrected data cube.
"""

import sys, os, shutil, obs
sys.path.append("..")
import sir
import definitions as d
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
from os.path import exists
from astropy.convolution import convolve_fft, convolve
from shutil import which

#############
# Help page #
#############
def help():
	print("correction_spectral_veil.py - Imports FTS data, shifts it, convolve it and compares it to the data")
	print("Usage: python correction_spectral_veil.py [OPTION]")
	print()
	sir.option("[1. Pos.]","config file")
	sys.exit()

#################################################################################################3

def argmin(x):
	"""
	Find the argument of the minimum in an multi-dimensional
	array. It seems to be faster than any numpy-built-in 
	function as shown in https://stackoverflow.com/questions/30180241/numpy-get-the-column-and-row-index-of-the-minimum-value-of-a-2d-array
	"""
	k = x.argmin()
	ncol = x.shape[1]
	return int(k/ncol), k%ncol

#################################################################################################3

def chi2(y_fit, y_obs):
     r"""
     Computes the merit-function $\chi^2$.
     
     Parameters
     ----------
     y_obs : numpy array
          Array with the observed values
     y_fit : numpy array
          Array with computed values, e.g. from a fit

     Return
     ------
     out : float
          $\chi^2$-value
     
     """
     return np.sum((y_fit-y_obs)**2)

#################################################################################################3

def gaussian(x, mean = 0, sigma = 1, norm = True):
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

	Return
	------
	out: float
		Value in the Gaussian function at position x
	
	"""
	g = np.exp(-(x - mean)**2 / (2 * sigma**2))
	if norm:
		g = g / (sigma * np.sqrt(2 * np.pi)) # Normalised so that int g = 1
	return g

#################################################################################################3

def lambda_0(ll,I):
	"""
	Determines the shifted lambda_0 due to convective blueshift

	Parameters
	----------
	ll : array
		Wavelengths in a given range
	I  : numpy.ndarray
		Stokes I in the given range

	Returns
	-------
	out : int	
		Index of the mean shifted spectral core
	out : float
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

def optimise_chi(nu, sigma, I, I_obs):
	"""
	Optimises two parameters by computing the chi2 merit function and finding
	the minima of the two parameters. The steps are the following
	 - Compute the chi2 values
	 - Compute the discrete minima in each parameter
	 - Perform two polynomial fits to find the optimal minimum
	 - Compute the uncertainty of the fit
	
	Parameter
	---------
	nu : array
		Array containing the nus which were used
	sigma : array
		Array containing the sigmas which were used
	I : ndarray
		Multidimensional array with the following format: I[nu, sigma, wavelength]. This are the fitted/computed data to be compared to the observed one
	I_obs : array
		Array containing the observed intensity
	
	Return
	------
	sigma_min : float
		Optimised sigma
	usigma_min : float
		Uncertainty of the optimised sigma
	nu_min : float
		Optimised nu
	unu_min : float
		Uncertainty of the optimised nu
	chis : ndarray
		Array containing all the chi2 for all possible combinations of nu and sigma
	
	"""
	chis = np.zeros(shape=(len(nu), len(sigma)))
	for n in range(len(nu)):
		for s in range(len(sigma)):
			chis[n,s] = np.log10(np.sum((I[n,s,:]-I_obs)**2))

	# Find the simple minima
	chi_min = argmin(chis)
	nu_min = nu[chi_min[0]]
	sigma_min = sigma[chi_min[1]]

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

def vac_to_air(wavelength):
	"""
	Computes the wavelength from vacuum to air by using (Ciddor, 1996) eq. 1.

	Note that this equation is only valid for wavelengths between 2300A and 16900A and
	it should be better in the infrared than Edlen.
	Note that the CO2 concentration is here assumed to be 450 ppm but it is most likely
	not constant in general.
	
	Parameters
	----------
	wavelength : float
		Wavelenght in Angstrom to be converted

	Return
	------
	out : float
		Corresponding wavelength in air in Angstrom
	"""
	sigma2 = (1e4 / wavelength)**2 # squared wavenumber in mum^-1
	refr    = (1 + 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2))
	return wavelength / refr

#################################################################################################3

def veil_correction(I_obs, nu, Ic = 1.):
	"""
	Correct the spectrum for the spectral veil by simly inverting
	the convoluted equation.

	Parameter
	---------
	I_obs : array
		Observed data
	nu : float
		Optimized fraction of spectral veil
	Ic : float, optional
		Continuum intensity of the FTS. Default is "1.0".
	"""
	
	return (I_obs - nu*Ic) / (1 - nu)

#################################################################################################3

def correct_spectral_veil(conf, stokes = ''):
	"""
	Correct spectral veil

	Parameter
	---------
	config : dict
		Dictionary with all the information from the config file
	stokes : array, optional
		Normalised data from the data cube ('' means to load the data from the config). Default: ''

	"""
	# Do not do anything
	if conf['fts_file'] == '':
		print("[STATUS] No spectral veil correction")
		if not exists(os.path.join(conf["path"], conf["cube_inv"])):
			if exists(os.path.join(conf["path"], conf["cube"]) + d.end_norm):
				print("[STATUS] Copying the created normalised data cube to the file used for the inversion.")
				shutil.copy(os.path.join(conf['path'],conf['cube']) + d.end_norm, os.path.join(conf['path'],conf['cube_inv'])) 
			elif exists(os.path.join(conf["path"], conf["cube"])):
				print("[STATUS] Copying the non-preprocessed data cube to the file used for the inversion.")
				shutil.copy(os.path.join(conf['path'],conf['cube']), os.path.join(conf['path'],conf['cube_inv'])) 
			else:
				print(f"[ERROR] File {conf['cube_inv']} does not exist and is not created. Are there data in the selected path?")
		return
	print("[STATUS] Correct spectral veil ...")
	if exists(os.path.join(conf['path'],conf['cube_inv'])):
		temp = ''
		print("[WARN] The data cube used for the inversion already exists and spectral veil correction is selected.")
		while temp != 'y' and temp != 'n':
			temp = input("Do you want to overwrite it and continue? [y/n] ")
		if temp == 'n':
			print("Abort (Consider changing the config file)")
			sys.exit(1)
	##########################
	#	Plot settings		#
	########################## 
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/mml.mplstyle'):
		plt.style.use(dirname + '/mml.mplstyle')
		# if dvipng is not installed, dont use latex
		if which('dvipng') is None:
			plt.rcParams["text.usetex"] = "False"
			plt.rcParams["font.family"] = 'sans-serif'
			plt.rcParams["mathtext.fontset"] = 'dejavuserif'
	elif "mml" in plt.style.available:
		plt.style.use('mml')
		# if dvipng is not installed, dont use latex
		if which('dvipng') is None:
			plt.rcParams["text.usetex"] = "False"
			plt.rcParams["font.family"] = 'sans-serif'
			plt.rcParams["mathtext.fontset"] = 'dejavuserif'
	else:
		plt.rcParams["savefig.format"] = "pdf"

	######################################################################################
	#						    DEFINE VARIABLES						    #
	######################################################################################
	b = 1.5 # For the used range in the spectrum
	sigma = np.arange(10,150,1) * 1e-3
	nu    = np.arange(0,31,1) * 1e-2

	######################################################################################
	#						    INITIALIZATION							 #
	######################################################################################
	# Load data
	if stokes == '':
		if len(conf['quiet_sun']) > 1:
			print("-------> Load normalised data ...")
			stokes = obs.load_data(conf, add = d.end_norm)
		else:
			print("[STATUS] Load data ...")
			stokes = obs.load_data(conf)
			
	filename_fts  = conf['fts_file']
	path = conf['path']

	# Read quiet sun file
	x1	  = conf['quiet_sun'][0]		# Lower limit for region in x
	x2	  = conf['quiet_sun'][1]+1	# Upper limit for region in x
	y1	  = conf['quiet_sun'][2]		# Lower limit for region in y
	y2	  = conf['quiet_sun'][3]+1	# Upper limit for region in y

	# FTS
	data = np.loadtxt(filename_fts)
	wn = np.flip(data[:,0]) # Wavenumber in cm^-1
	i1 = np.flip(data[:,1]) # Intensity corrected


	# GRIS
	ll_gris = np.load(os.path.join(conf['path'],conf['waves']))

	######################################################################################
	#						   SHIFT OF SPECTRUM						    #
	######################################################################################
	print('-------> Shift spectrum ...')
	###########
	#    FTS  #
	###########
	# Literature wavelength in air where the line core is expected in air
	if conf['instrument'] in d.ll_lit:
		ll_lit = d.ll_lit[conf['instrument']]
	else:
		ll_lit = float(input("Used instrument not defined in definitions. Type the literature wavelength in air, used for the spectral veil correction: "))
	#ll_lit = 15648.515 # From A new Multiplet table for Fe I; imported from definitions

	# Compute wavelength in vacuum and in air
	ll = 1e8 / wn # in wavelength in A
	ll_ciddor = vac_to_air(ll) # Ciddor, 1996 method

	#################
	#    Data		 #
	#################
	# Compute average intensity in quiet sun in data and shift it to lambda_0
	i_gris = np.mean(stokes[x1:x2,y1:y2,0,:], axis=(0,1))

	# Correct convective blueshift / Shift in GRIS data
	ll1_lit = np.argmin(abs(ll_gris-ll_lit))		 # Find position at ll_lit defined in definitions.py
	border = 40								 # Find minima around this value
	ll1, l0 = lambda_0(ll_gris[ll1_lit-border:ll1_lit+border], stokes[x1:x2,y1:y2,0,ll1_lit-border:ll1_lit+border]) # Determine the minima position and value of the shifted GRIS spectrum
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
	i_conv	= i1		 [np.argmin(abs(ll_ciddor + b)) : np.argmin(abs(ll_ciddor - b))]

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
	np.save(os.path.join(conf['path'],d.veil_parameters), [nu_min, sigma_min]) # sigma in Angstrom

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
	c = ax.contourf(sigma*1e3, nu*100, chis, levels = np.linspace(vmin, vmax, 40), cmap = 'ocean', vmin = vmin, vmax = vmax)
	ax.set_ylabel(r"$\nu$ $[\%]$")
	ax.set_xlabel(r"$\sigma$ [m\AA]")
	cbar = fig.colorbar(c,ticks=np.arange(vmin, vmax+0.2, 0.2), location = 'top')
	cbar.set_label(label=r"$\log(\chi^2)$", loc = 'center', labelpad = 15)

	plt.text(120, 25, r"$\hat{\sigma} =$ " + "(%.1f ± %.1f)"
					% (round(sigma_min/1e-3,1), round(usigma_min / 1e-3 + 0.05,1)) + r" m\AA"
				 + '\n' +
				 r"$\hat{\nu} =$ " + "(%.1f ± %.1f) "
					% (round(nu_min/1e-2,1), round(unu_min / 1e-2 + 0.05,1))  + r"$\%$"
				 + '\n' +
				 r"$\hat{\chi}^2 =$ " + "%.3f"
					% (chi_min), fontsize=12)
	plt.savefig(os.path.join(path,"veil_parameter_space"))


	######################################################################################
	#					    Plot the best convolved intensity				  #
	######################################################################################
	fig, ax = plt.subplots()
	ax.plot(ll_conv, i_conv, "x", label = r'$I^{FTS}$')
	if conf['instrument'] == 'GRIS':
		ax.plot(ll_conv_g, I_obs, '+', label = r'$I_{\mathrm{qs}}^{\mathrm{GRIS}}$')
	else:
		ax.plot(ll_conv_g, I_obs, '+', label = r'$I_{\mathrm{qs}}^{\mathrm{data}}$')
	ax.plot(LL_conv, I_conv_best, label = r"$I^{\mathrm{fit}}$")

	ax.set_xlabel(r"$\Delta \lambda^{\mathrm{air}}$ [\AA]")
	ax.set_ylabel(r"Intensity $I/I_c$")
	ax.set_xlim(-b,b)

	ax.legend(fontsize=15)
	plt.savefig(os.path.join(path, "veil_intensity"))

	######################################################################################
	#					START OF THE ACTUAL CORRECTION					    #
	######################################################################################
	print('[STATUS] Correct data ...')

	nx, ny, ns, nw = stokes.shape[0], stokes.shape[1], stokes.shape[2], stokes.shape[3]
	data = np.empty((nx, ny, ns, nw))

	# Correct for spectral veil
	data[:,:,0,:] = veil_correction(stokes[:,:,0,:], nu_min)
	data[:,:,1,:] = stokes[:,:,1,:] # No correction needed in Q
	data[:,:,2,:] = stokes[:,:,2,:] # No correction needed in U
	data[:,:,3,:] = stokes[:,:,3,:] # No correction needed in V

	# To save data convert to float32
	data = data.astype(np.float32)

	print("-------> Saving data (this might take a while) ...")
	if ".npy" in conf['cube_inv']:
		np.save(os.path.join(conf['path'],conf['cube_inv']), data)
	else:
		# Write the merged data cube
		example = fits.open(os.path.join(path,conf['cube']))
		header = example[0].header
		hdu = fits.PrimaryHDU(data)
		hdu.header = header
		hdu.writeto(os.path.join(path,conf['cube_inv']),overwrite=True)

if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])	
	correct_spectral_veil(conf, stokes = '')






