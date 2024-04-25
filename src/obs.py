"""
Functions related to observations and to conversions between SIR format and numpy files
"""
import numpy as np
import sys, os, sir
from astropy.io import fits 
import definitions as d
from os.path import exists

##############################################################################################
def load_data(conf, filename = '', add = '', add_path = True):
	"""
	Load the data cube as given in the config
	
	Parameter
	---------
	config : dict
		Dictionary with all the configurations
	filename : string, optional
		If empty, filename from config is used. Default = ''
	add : string, optional
		Add this string before the extension to open corrected or
		normalised data cube. Default: ''
	add_path : bool, optional
		Add the path in the config. Default: True

	Return
	------
	numpy array

	"""
	if filename == '':
		filename = conf['cube']

	# Open norm or corrected data 
	filename = filename.replace('.fits',add + '.fits').replace('.npy',add + '.npy')

	if add_path:
		cube = os.path.join(conf['path'],filename)
	else:
		cube = filename

	if not exists(cube):
		print(f"[load data] File {cube} does not exist!")
		if add == d.end_norm:
			print(f"            Is the normalised data cube created?")
		sys.exit()
		
	if ".npy" in cube:
		data = np.load(cube)
	else:
		example = fits.open(cube)
		data = example[0].data
		data = data

	return data

##############################################################################################

def read_profile(profile, grid, line_file, waves):

	# Read grid and line file to transform the relative llambdas
	grid = sir.read_grid(grid)
	Lines = grid['Line']
	Mins  = grid['min']
	Steps = grid['step']
	Maxs  = grid['max']

	#Line file
	line = sir.read_line(line_file)
	wave = line['wavelength']
	line = line['Line']

	
	# Read the profile
	lls = np.array([])
	Is = np.array([])
	Qs = np.array([])
	Us = np.array([])
	Vs = np.array([])
	for i in Lines:
		ll, I, Q, U, V = sir.read_profile(profile, num = int(i[0]))
		lls = np.append(lls, np.array(ll) + wave[np.where(line == int(i[0]))]) # in abs. wavelength
		Is = np.append(Is, I)
		Qs = np.append(Qs, Q)
		Us = np.append(Us,U)
		Vs = np.append(Vs, V)

	return lls, Is, Qs, Us, Vs


##############################################################################################

def write_psf(conf, filename):
	"""
	Writes the spectral point spread function with the value from the spectral veil correction

	Parameter
	---------
	config : dict
		Dictionary containing all the information of the config file
	filename : string
		Filename under which the file is saved

	Return
	------
	None

	"""
	waves = np.load(os.path.join(conf['path'],conf['waves'])) # Wavelengths
	Delta_ll = waves[1]-waves[0]

	## Wavelength range to be printed in mA
	ll = np.arange(-(Delta_ll*1e3*15),(Delta_ll*1e3*16),Delta_ll*1e3)

	# Check if veil parameters are saved in the path
	if not exists(os.path.join(conf['path'],d.veil_parameters)):
		print(f"[write_psf] ERROR: Veil Parameters file {d.veil_parameters} does not exist. Was it corrected?")
		print(f"            No PSF file is created. Consider removing the parameter in the config file or provide the file.")
		sys.exit()

	# Load sigma in mA from the spectral veil correction
	sigma = np.load(os.path.join(conf['path'],d.veil_parameters))[0]*1e3


	# Compute Gaussian
	g = np.exp(-(ll - 0)**2 / (2 * sigma**2))
	g = g / (sigma * np.sqrt(2 * np.pi)) # Normalised so that int g = 1

	f = open(filename, 'w')
	for num1,num2 in zip(ll,g):
		f.write(f" {num1:>9.4f}   {num2:>10.4E}\n")







