"""
Library for repeating functions for analyzing and/or plotting SIR data
"""
import numpy as np 
import matplotlib.pyplot as plt 
import fileinput, time, sys, os, definitions as d
from os.path import exists
import obs

######################################################################
def angstrom_to_pixel(waves, range_wave):
	"""
	Converts angstrom to pixel for a given waves file.
	
	Parameter
	---------
	waves : numpy array
		Array containing the real wavelengths in Angstrom
	range_wave : numpy array
		Array containing the ranges in Angstrom

	Return
	------
	Array with the ranges in pixel for the given waves array
	"""
	range_wave_new = np.copy(range_wave) # Otherwise the config is overwritten
	for i in range(len(range_wave)):
		if range_wave[i][0] >= waves[0]:
			range_wave_new[i][0] = np.abs(waves-range_wave[i][0]).argmin()
			range_wave_new[i][1] = np.abs(waves-range_wave[i][1]).argmin()
	range_wave_new = range_wave_new.astype(int) # Convert to integer

	return range_wave_new

def pixel_to_angstrom(waves, range_wave):
	"""
	Converts pixel to angstrom for a given waves file.
	
	Parameter
	---------
	waves : numpy array
		Array containing the real wavelengths in Angstrom
	range_wave : numpy array
		Array containing the ranges in Angstrom

	Return
	------
	Array with the ranges in Angstrom for the given waves array
	"""
	range_wave_new = np.copy(range_wave) # Otherwise the config is overwritten
	for i in range(len(range_wave)):
		if range_wave[i][0]  < waves[0]:
			range_wave_new[i][0] = waves[range_wave[i][0]]
			range_wave_new[i][1] = waves[range_wave[i][1]]

	return range_wave_new



######################################################################

def read_config(filename, check = True, change_config = False):
	"""
	Reads a config file for the inversion
	
	Parameter
	---------
	filename : string
		Path of the control file
	check : bool, optional
		Check if file exists (Default: True)
	change_config : bool, optional
		config file is read to be changed (=> Do not try to load anything) (Default: False)
	

	Return
	------
	Dict : dict
		Dict containing all the information from the config file
	
	"""
	if not exists(filename):
		print("[read_config] " + filename + " does not exist.")
		sys.exit(1)

	# Load data
	data = np.genfromtxt(filename, delimiter=':', comments='#', dtype=str)

	# Remove spaces in first column
	for i in data:
		# Check if the value is empty
		if len(i[0].replace(' ','')) == 0:
			i[0] = ''
			continue
		while i[0][0] == ' ':
			i[0] = i[0][1:]
		while i[0][-1] == ' ':
			i[0] = i[0][:-1]

	# Remove spaces in second column
	for i in data:
		# Check if the value is empty
		if len(i[1].replace(' ','')) == 0:
			i[1] = ''
			continue
		while i[1][0] == ' ':
			i[1] = i[1][1:]
		while i[1][-1] == ' ':
			i[1] = i[1][:-1]

		i[0] = i[0].replace('  ','')   # replace double spaces

	# Get the name of the file
	if "/" in filename:
		filename = filename[filename.rfind('/')+1:]

	# Create dictionary
	Dict = {
		"filename" : filename,
		"vmacro" : "0.1000",
	}

	for i in data:
		if i[0] != '':
			Dict[i[0]] = i[1]

	if Dict["inv_out"] == "":
		print("[read_config] No value in inv_out")
	if Dict["line"] == "":
		print("[read_config] No value in line")
	if Dict["atoms"] == "":
		print("[read_config] No value in atoms")
	if Dict["cycles"] == '':
		print("[read_config] No value in cycles")
	if Dict["abundance"] == '':
		print("[read_config] No value in abundance")

	# Transform the information into lists or in different types than string
	if Dict['mode'] == "1C" or Dict['mode'] == "2C":
		Dict['map'] = np.array([int(i) for i in Dict["map"].split(',')], dtype=int)
		Dict["quiet_sun"] = np.array([int(i) for i in Dict["quiet_sun"].split(',')])
		
	if Dict['mode'] == "MC":
		Dict['num'] = int(Dict['num'])

	# Convert the ranges into integers or floats
	temp = [i. split(',') for i in Dict["range_wave"].split(';')]
	for i in range(0,len(temp)):
		temp[i][0] = float(temp[i][0])
		temp[i][1] = float(temp[i][1])
	
	Dict["range_wave"] = np.array(temp)

	Dict['atoms'] = Dict["atoms"].split(';')	# Atoms

	Dict["random_guess"] = int(Dict["random_guess"])

	Dict["random_pars"] = Dict["random_pars"].replace(" ",'').split(",")
	
	Dict["cycles"] = int(Dict["cycles"])

	Dict["weights"] = Dict["weights"].split(',')

	# Check if range_wave fits the atoms
	if len(Dict['atoms']) != len(Dict['range_wave']):
		print("[read_config] The number of lines in 'atoms' do not fit the given ranges in 'range_wave'! Abort...")
		sys.exit()

	# Check if some files exists
	if check:
		if not exists(Dict['path']):
			print(f"[read_config] {Dict['path']} does not exist.")
		if not exists(os.path.join(Dict['path'],Dict['model'])):
			print(f"[read_config] {Dict['model']} does not exist.")
		if Dict['mode'] == "1C" or Dict['mode'] == "2C":
			if not exists(os.path.join(Dict['path'],Dict['fts_file'])) and Dict['preprocess'] == "1":
				print(f"[read_config] {Dict['fts_file']} does not exist.")
	
	# If Path is relative, change to absolute
	if Dict['path'][0:2] == "./":
		if Dict['path'] == "./":
			Dict['path'] = os.getcwd()
		else:
			Dict['path'] = os.getcwd() + Dict[2:]
	
	# Correction for old version
	if "lim_azimuth" in Dict:
		Dict["lim_phi"] = Dict["lim_azimuth"]
	if "lim_azimuth1" in Dict:
		Dict["lim_phi1"] = Dict["lim_azimuth1"]
	if "lim_azimuth2" in Dict:
		Dict["lim_phi2"] = Dict["lim_azimuth2"]

	return Dict

######################################################################
def read_control(filename):
	"""
	Reads a control file in the scheme SIR expects it.
	
	Parameter
	---------
	filename : string
		Path of the control file
	
	Return
	------
	Dict : dict
		Dict containing all the information from the control file
	
	"""
	if not exists(filename):
		print("[read_control] " + filename + " does not exist.")
		sys.exit(1)

	# Load data
	data = np.genfromtxt(filename, delimiter=':', comments='!', dtype=str)
	for i in range(len(data)):
		data[i][0] = data[i][0].replace('(*)','') # replace the (*)
		data[i][0] = data[i][0].replace('  ','') # replace double spaces
	
	Dict = dict()
	for i in data:
		Dict[i[0]] = i[1]
	return Dict

######################################################################
def read_model(filename):
	"""
	Reads a model file and returns all parameters
	
	Parameter
	---------
	filename : string
		String containing the path of the file

	Return
	------
	log_tau : numpy.array
		Log tau
	T : numpy.array
		Temperature in K
	Pe : numpy.array
		Electron pressure in dyn/cm^2
	v_micro : numpy.array
		Microturbulence velocity in cm/s
	B : numpy.array
		Magnetic field strength in Gauss
	vlos : numpy.array
		Line-of-sight velocity in cm/s
	inc : numpy.array
		Inclination in deg
	azimuth : numpy.array
		Azimuth angle in deg
	z : numpy.array, optional
		Height in km
	Pg : numpy.array, optional
		Gas pressure in dyn/cm^2
	rho : numpy.array, optional
		Density in g/cm^3
	"""	
	# Open file

	data = np.loadtxt(filename, skiprows=1).transpose()

	# Store data
	log_tau	= data[0]
	T		= data[1]
	Pe	 	= data[2]
	v_micro 	= data[3]
	B	 	= data[4]
	vlos     	= data[5]
	inc  	= data[6]
	azimuth	= data[7]
	z	 	= np.zeros(shape=data[7].shape)
	Pg	 	= np.zeros(shape=data[7].shape)
	rho		= np.zeros(shape=data[7].shape)

	# If data contains more rows, add the data
	if len(data) > 8:
		z   = data[8]
		Pg  = data[9]
		rho = data[10]
	
	return log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho

	
######################################################################
def read_profile(filename, num = 0):
	"""
	Reads the first LINE data from a profile computed by SIR
	
	Parameter
	---------
	filename : string
		String containing the path of the file
	num : int, optional
		Number of the line which is loaded. Default: 0 (use first one from line)

	Return
	------
	ll : numpy.array
		Wavelengths in A
	I : numpy.array
		Stokes I
	Q : numpy.array
		Stokes Q
	U : numpy.array
		Stokes U
	V : numpy.array 
		Stokes V
	"""
	num = int(num) # num must be an integer
	
	data = np.loadtxt(filename).transpose()
	line = data[0].astype(int)
	ll = data[1].astype(np.float64)
	I  = data[2].astype(np.float64)
	Q  = data[3].astype(np.float64)
	U  = data[4].astype(np.float64)
	V  = data[5].astype(np.float64)

	# Only return values for line x
	if num == 0:
		num = line[0]

	if num != -1:
		ll = ll[line == num]
		I  = I [line == num]
		Q  = Q [line == num]
		U  = U [line == num]
		V  = V [line == num]
		return ll/1000, I, Q, U, V
	
	else:
		return np.array(line), np.array(ll/1000), np.array(I), np.array(Q), np.array(U), np.array(V)

######################################################################
def read_grid(filename):
	"""
	Reads the grid file
	
	Parameter
	---------
	filename : string
		File to be read
	
	Return
	-------
	dict : Dictionary
		Dict. with 'Line', 'min', 'step' and 'max' in it
	"""
	# Open the file and read lines
	with open(filename) as f:
		strings = f.readlines()

	# Remove last line if it contains no information
	if ":" not in strings[-1]:
		strings = strings[:-1]

	# Remove leading spaces
	for i in range(len(strings)):
		while(strings[i][0] == ' '):
			strings[i] = strings[i][1:]
	
	# Create arrays
	Line	  = []
	ll_min	= np.empty(len(strings))
	ll_step    = np.empty(len(strings))
	ll_max	= np.empty(len(strings))

	
	# Remove multiple spaces and split 
	strings = [i.replace(':', ': ')   for i in strings] # Add space to make sure it is split correctly
	while any('  ' in x for x in strings):
		strings = [i.replace('  ', ' ')   for i in strings]
	strings = [i.split(' ') for i in strings]

	# Fill arrays with data	
	for i in range(len(strings)):
		Line.append(strings[i][0][:-1].split(','))
		ll_min[i]  = float(strings[i][1].replace(",",""))
		ll_step[i] = float(strings[i][2].replace(",",""))
		ll_max[i]  = float(strings[i][3].replace(",",""))

	Dict = {
			'Line'		: Line,
			'min'		: ll_min,
			'step'		: ll_step,
			'max'		: ll_max,
		  }
	
	return Dict

######################################################################
def read_line(filename):
	"""
	Reads the line file

	Parameter
	---------
	filename : string
		File to be read
	
	Return
	-------
	dict : Dictionary
		Dict. with 'Line', 'Ion', 'wavelength', 'factor', 'Exc_Pot', log_gf',
		'Transition', 'alpha' and 'sigma' in it
	"""
	# Open the file and read lines
	with open(filename) as f:
		strings = f.readlines()

	# Remove last line if it contains no information
	if "=" not in strings[-1]:
		strings = strings[:-1]

	# Remove leading spaces
	for i in range(len(strings)):
		while(strings[i][0] == ' '):
			strings[i] = strings[i][1:]
	
	# Create arrays
	Line	   = np.empty(len(strings), dtype = int)
	Ion	    = ["" for i in strings]
	ll		= np.empty(len(strings))
	factor	 = np.empty(len(strings))
	Exc_Pot	= np.empty(len(strings))
	log_gf	 = np.empty(len(strings))
	Transition  = np.empty(len(strings), dtype=str)
	alpha	  = np.empty(len(strings))
	sigma	  = np.empty(len(strings))
	
	# Remove multiple spaces
	while any('  ' in x for x in strings):
		strings = [i.replace('  ', ' ')   for i in strings]

	# Fill arrays with data	
	strings = [i.split(' ') for i in strings]
	for i in range(len(strings)):
		split = strings[i][0].split('=')
		Line[i]    = int(split[0])
		Ion[i]	= split[1] + ' ' + strings[i][1]
		ll[i]	= strings[i][2]
		factor[i]  = strings[i][3]
		Exc_Pot[i] = strings[i][4]
		log_gf[i]  = strings[i][5]
		Transition[i]  = strings[i][6] + strings[i][7] + strings[i][8] + strings[i][9]
		alpha[i]   = strings[i][10]
		sigma[i]   = strings[i][11]

	Dict = {
			'Line'		: Line,
			'Ion'		: Ion,
			'wavelength'	: ll,
			'factor'		: factor,
			'Exc_Pot'		: Exc_Pot,
			'log_gf'		: log_gf,
			'Transition'	: Transition,
			'alpha'	 : alpha,
			'sigma'	 : sigma
		  }
	
	return Dict

######################################################################
def list_to_string(temp, let = ','):
	"""
	Convert a list to a string

	Parameter
	---------
	temp : list
	let : str
		Letter which is added as a separation

	Return
	------
	string with the information from the list
	"""
	temp1 = ''
	for i in range(len(temp)):
		temp1 += str(temp[i])
		if i + 1 < len(temp):
				temp1 += let
	return temp1

################################################################################

def write_config_1c(File, conf):
	"""
	Writes a config file with the information provided as a dictionary for the mode 1C

	Parameters
	----------
	File : string
		Save path
	conf : dict
		Dictionary with all the informations
	"""
	print("[write_config] Manually added comments will be overwritten? 1s to abort left ...")
	time.sleep(2)
	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	for i in range(len(temp)):
		for j in range(len(temp[i])):
			range_wave += str(temp[i][j])
			if j +1 < len(temp[i]):
				range_wave += ','
		if i < len(temp)-1:
			range_wave += ';'

	atoms = list_to_string(conf["atoms"], ";")
	Map = list_to_string(conf["map"])
	weights = list_to_string(conf["weights"])
	quiet_sun = list_to_string(conf["quiet_sun"])
	random_pars = list_to_string(conf["random_pars"])
	

	with open(File, 'w') as f:
		f.write("# This is the config file, generated with create_config.py\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")

		f.write(f"# \n")
		f.write(f"# Stuff from the data\n")
		f.write(f"# \n")
		f.write(f"cube : {conf['cube']} # Data cube name (bin) used for preprocessing data if 'preprocess' is 1\n")
		f.write(f"cube_inv : {conf['cube_inv']} # Data cube name used for the inversion (bin)\n")
		f.write(f"map : {Map} # Pixels to be considered as a list\n")
		f.write(f"instrument : {conf['instrument']} # Instrument used (GRIS, Hinode or empty)\n")
		f.write(f"shift_wave : {conf['shift_wave']} # Shift the wavelength grid when waves file is created in mA\n")

		f.write("#\n")
		f.write("# Data Preprocessing when main.py is executed\n")
		f.write("#\n")
		f.write(f"preprocess : {conf['preprocess']} # Preprocess data (1 = True, 0 = False)\n")
		f.write(f"quiet_sun : {quiet_sun} # Quiet sun region for normalization as a list (0 => already normalised)\n")
		f.write(f"fts_file : {conf['fts_file']} # Absolute fts file, blank = do not correct spectral veil\n")

		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"model : {conf['model']} # Base Model for guess\n")
		f.write(f"range_wave : {range_wave} # Ranges of wavelengths (pixel or Angstrom) to be considered min1,max1;min2,max2;... First pair belongs to first line in Grid file, etc.\n")
		f.write(f"inv_out : {conf['inv_out']} # Prefix of output of the inversion files\n")
		f.write(f"chi2 : {conf['chi2']} # Output of the chi2 values (npy)\n")
		f.write(f"line : {conf['line']} # Line file\n")
		f.write(f"atoms : {atoms} # Atoms used, ; defines a new line\n")
		f.write(f"guess : {conf['guess']} # Use a bin file as initial guesses, blank use base model\n")
		f.write(f"psf : {conf['psf']} # .dat file (if it does not exist, computed from spectral veil parameter), blank=not used\n")

		f.write(f"# \n")
		f.write(f"# Control file\n")
		f.write(f"# \n")
		f.write(f"cycles : {conf['cycles']} # Number of cycles\n")
		f.write(f"weights : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp : {conf['nodes_temp']} # Nodes in T\n")
		f.write(f"nodes_magn : {conf['nodes_magn']} # Nodes in B\n")
		f.write(f"nodes_vlos : {conf['nodes_vlos']} # Nodes in vlos\n")
		f.write(f"nodes_gamma : {conf['nodes_gamma']} # Nodes in gamma\n")
		f.write(f"nodes_phi : {conf['nodes_phi']} # Nodes in phi\n")
		f.write(f"vmacro : {conf['vmacro']} # Macroturbulence velocity\n")
		f.write(f"mu_cos : {conf['mu_cos']} # mu = cos theta\n")
		f.write(f"abundance : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition\n")

		f.write(f"# \n")
		f.write(f"# Radomisation Settings\n")
		f.write(f"# \n")
		f.write(f"random_guess : {conf['random_guess']} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars : {random_pars} # Randomise these parameters in the file(s) below\n")
		f.write(f"lim_B : {conf['lim_B']} # Limits for the randomisation in B in G\n")
		f.write(f"lim_vlos : {conf['lim_vlos']} # Limits for the randomisation in vlos in cm/s\n")
		f.write(f"lim_gamma : {conf['lim_gamma']} # Limits for the randomisation in the inclination in deg\n")
		f.write(f"lim_phi : {conf['lim_phi']} # Limits for the randomisation in the azimuth in deg")

######################################################################

def write_config_2c(File, conf):
	"""
	Writes a config file with the information provided as a dictionary

	Parameters
	----------
	File : string
		Save path
	conf : dict
		Dictionary with all the informations
	"""
	print("[write_config] Manually added comments will be overwritten? 1s to abort left ...")
	time.sleep(2)

	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	for i in range(len(temp)):
		for j in range(len(temp[i])):
			range_wave += str(temp[i][j])
			if j +1 < len(temp[i]):
				range_wave += ','
		if i < len(temp)-1:
			range_wave += ';'

	atoms = list_to_string(conf["atoms"], ";")
	Map = list_to_string(conf["map"])
	weights = list_to_string(conf["weights"])
	quiet_sun = list_to_string(conf["quiet_sun"])
	random_pars = list_to_string(conf["random_pars"])
	

	with open(File, 'w') as f:
		f.write("# This is the config file, generated with create_config.py\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")

		f.write(f"# \n")
		f.write(f"# Stuff from the data\n")
		f.write(f"# \n")
		f.write(f"cube : {conf['cube']} # Data cube name (bin) used for preprocessing data if 'preprocess' is 1\n")
		f.write(f"cube_inv : {conf['cube_inv']} # Data cube name used for the inversion (bin)\n")
		f.write(f"map : {Map} # Pixels to be considered as a list (0 means all pixels)\n")
		f.write(f"instrument : {conf['instrument']} # Instrument used (GRIS, Hinode or empty)\n")
		f.write(f"shift_wave : {conf['shift_wave']} # Shift the wavelength grid when waves file is created in mA\n")

		f.write("#\n")
		f.write("# Data Preprocessing when main.py is executed\n")
		f.write("#\n")
		f.write(f"preprocess : {conf['preprocess']} # Preprocess data (1 = True, 0 = False)\n")
		f.write(f"quiet_sun : {quiet_sun} # Quiet sun region for normalization as a list (0 => already normalised)\n")
		f.write(f"fts_file : {conf['fts_file']} # Absolute fts file, blank = do not correct spectral veil\n")


		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"model1 : {conf['model1']} # Base Model 1 for guess\n")
		f.write(f"model2 : {conf['model2']} # Base Model 2 for guess\n")
		f.write(f"range_wave : {range_wave} # Ranges of wavelengths (pixel or Angstrom) to be considered min1,max1;min2,max2;... First pair belongs to first line in Grid file, etc.\n")
		f.write(f"inv_out : {conf['inv_out']} # Prefix of output of the inversion files\n")
		f.write(f"chi2 : {conf['chi2']} # Output of the chi2 values (npy)\n")
		f.write(f"line : {conf['line']} # Line file\n")
		f.write(f"atoms : {atoms} # Atoms used, ; defines a new line\n")
		f.write(f"guess1 : {conf['guess1']} # Use a bin file as initial guesses, blank use base model 1\n")
		f.write(f"guess2 : {conf['guess2']} # Use a bin file as initial guesses, blank use base model 2\n")
		f.write(f"psf : {conf['psf']} # .dat file (if it does not exist, computed from spectral veil parameter), blank=not used\n")

		f.write(f"# \n")
		f.write(f"# Control file\n")
		f.write(f"# \n")
		f.write(f"cycles : {conf['cycles']} # Number of cycles\n")
		f.write(f"weights : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp1 : {conf['nodes_temp1']} # Nodes 1 in T\n")
		f.write(f"nodes_magn1 : {conf['nodes_magn1']} # Nodes 1 in B\n")
		f.write(f"nodes_vlos1 : {conf['nodes_vlos1']} # Nodes 1 in vlos\n")
		f.write(f"nodes_gamma1 : {conf['nodes_gamma1']} # Nodes 1 in gamma\n")
		f.write(f"nodes_phi1 : {conf['nodes_phi1']} # Nodes in 1 phi\n")
		f.write(f"nodes_temp2 : {conf['nodes_temp2']} # Nodes 2 in T\n")
		f.write(f"nodes_magn2 : {conf['nodes_magn2']} # Nodes 2 in B\n")
		f.write(f"nodes_vlos2 : {conf['nodes_vlos2']} # Nodes 2 in vlos\n")
		f.write(f"nodes_gamma2 : {conf['nodes_gamma2']} # Nodes 2 in gamma\n")
		f.write(f"nodes_phi2 : {conf['nodes_phi2']} # Nodes 2 in phi\n")
		f.write(f"mu_cos : {conf['mu_cos']} # mu = cos theta\n")
		f.write(f"abundance : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition\n")

		f.write(f"# \n")
		f.write(f"# Radomisation Settings\n")
		f.write(f"# \n")
		f.write(f"random_guess : {conf['random_guess']} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars : {random_pars} # Randomise these parameters in the file(s) below\n")
		f.write(f"lim_B1 : {conf['lim_B1']} # Limits 1 for the randomisation in B in G\n")
		f.write(f"lim_vlos1 : {conf['lim_vlos1']} # Limits 1 for the randomisation in vlos in cm/s\n")
		f.write(f"lim_gamma1 : {conf['lim_gamma1']} # Limits 1 for the randomisation in the inclination in deg\n")
		f.write(f"lim_phi1 : {conf['lim_phi1']} # Limits 1 for the randomisation in the azimuth in deg\n")
		f.write(f"lim_B2 : {conf['lim_B2']} # Limits 1 for the randomisation in B in G\n")
		f.write(f"lim_vlos2 : {conf['lim_vlos2']} # Limits 2 for the randomisation in vlos in cm/s\n")
		f.write(f"lim_gamma2 : {conf['lim_gamma2']} # Limits 2 for the randomisation in the inclination in deg\n")
		f.write(f"lim_phi2 : {conf['lim_phi2']} # Limits 2 for the randomisation in the azimuth in deg")

######################################################################

def write_config_mc(File, conf):
	"""
	Writes a config file with the information provided as a dictionary

	Parameters
	----------
	File : string
		Save path
	conf : dict
		Dictionary with all the informations
	"""
	print("[write_config] Note that manually added comments will be overwritten! 1s left to abort ...")
	time.sleep(2)
	
	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	for i in range(len(temp)):
		for j in range(len(temp[i])):
			range_wave += str(temp[i][j])
			if j +1 < len(temp[i]):
				range_wave += ','
		if i < len(temp)-1:
			range_wave += ';'

	atoms = list_to_string(conf["atoms"], ";")
	weights = list_to_string(conf["weights"])
	random_pars = list_to_string(conf["random_pars"])
	
	with open(File, 'w') as f:
		f.write(f"# Config file\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")
		f.write(f"# \n")
		f.write(f"# General Stuff\n")
		f.write(f"# \n")
		f.write(f"num : {conf['num']} # Number of Models\n")
		f.write(f"instrument : {conf['instrument']} # Instrument used (GRIS, Hinode or empty)\n")
		f.write(f"model : {conf['model']} # Base Model for guess\n")
		f.write(f"atoms : {atoms} # Atoms to be used in Grid file\n")	
		f.write(f"range_wave : {range_wave} # Ranges of wavelengths in mA to be considered min1,step1,max1;min2,step2,max2;... First pair belongs to first line in Grid file, etc.\n")
		f.write(f"random_guess : {conf['random_guess']} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars : {random_pars} # Randomise these parameters for the guess\n")
		f.write(f"model_pars : {conf['model_pars']} # Randomise these parameters while creating models as a list\n")
		
		f.write(f"#\n")
		f.write(f"# Creating Models and Synthesis\n")
		f.write(f"#\n")
		f.write(f"model_nodes : {conf['model_nodes']} # Create models with 2 or 3 nodes\n")
		f.write(f"model_out : {conf['model_out']} # Output file of the created models\n")
		f.write(f"syn_out : {conf['syn_out']} # Output of the synthesis profiles and models\n")
		f.write(f"noise_out : {conf['noise_out']} # Output of the noise profiles\n")
		f.write(f"noise_I : {conf['noise_I']} # Noise in I\n")		
		f.write(f"noise_Q : {conf['noise_Q']} # Noise in Q\n")
		f.write(f"noise_U : {conf['noise_U']} # Noise in U\n")
		f.write(f"noise_V : {conf['noise_V']} # Noise in V\n")
		f.write(f"create_B : {conf['create_B']} # The limits for the first and last node in B\n")
		f.write(f"create_vlos : {conf['create_vlos']} # The limits for the first and last node in vlos\n")
		f.write(f"create_gamma : {conf['create_gamma']} # The limits for the first and last node in gamma\n")
		f.write(f"create_phi : {conf['create_phi']} # The limits for the first and last node in phi\n")
		f.write(f"create_points : {conf['create_points']} # At this log tau points the models are interpolated with splines (increasing), 2 or 3 values for 2 or 3 nodes\n")

		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"inv_out : {conf['inv_out']} # Prefix of the output of the inversion files\n")
		f.write(f"chi2 : {conf['chi2']} # Output of the chi2 values (npy)\n")
		f.write(f"line : {conf['line']} # Line file\n")
		f.write(f"guess : {conf['guess']} # Use a bin file as initial guesses, blank use base model\n")

		f.write(f"# \n")
		f.write(f"# Control file\n")
		f.write(f"# \n")
		f.write(f"cycles : {conf['cycles']} # Number of cycles\n")
		f.write(f"weights : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp : {conf['nodes_temp']} # Nodes in T\n")
		f.write(f"nodes_magn : {conf['nodes_magn']} # Nodes in B\n")
		f.write(f"nodes_vlos : {conf['nodes_vlos']} # Nodes in vlos\n")
		f.write(f"nodes_gamma : {conf['nodes_gamma']} # Nodes in gamma\n")
		f.write(f"nodes_phi : {conf['nodes_phi']} # Nodes in phi\n")
		f.write(f"vmacro : {conf['vmacro']} # Macroturbulence velocity\n")
		f.write(f"abundance : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition\n")

		f.write(f"# \n")
		f.write(f"# Radomisation Settings\n")
		f.write(f"# \n")
		f.write(f"lim_B : {conf['lim_B']} # Limits for the randomisation in B in G\n")
		f.write(f"lim_vlos : {conf['lim_vlos']} # Limits for the randomisation in vlos in cm/s\n")
		f.write(f"lim_gamma : {conf['lim_gamma']} # Limits for the randomisation in the inclination in deg\n")
		f.write(f"lim_phi : {conf['lim_phi']} # Limits for the randomisation in the azimuth in deg")

######################################################################

def write_config(File, conf):
	"""
	Writes a config file with the information provided as a dictionary for the mode 1C

	Parameters
	----------
	File : string
		Save path
	conf : dict
		Dictionary with all the informations
	"""
	if conf["mode"] == "MC":
		write_config_mc(File, conf)
	elif conf["mode"] == "1C":
		write_config_1c(File, conf)
	elif conf["mode"] == "2C":
		write_config_2c(File, conf)
	else:
		print("[write_config] Mode is not defined and config file cannot be written.")

######################################################################

def write_grid(conf, waves, filename = 'Grid.grid'):
	"""
	Writes the Grid file with data from the config file

	Parameter
	---------
	config : dict
		Dictionary containing all the information from the config file
	filename : string, optional
		String containing the name of the Grid file. Default: Grid.grid
	"""

	# Load data from config
	range_wave = conf['range_wave']
	line = read_line(os.path.join(conf['path'],conf['line']))
	atoms = conf["atoms"]

	# Change the ranges if it is given in angstroms, if needed
	range_wave = pixel_to_angstrom(waves,range_wave)

	# Define minimum, step and maximum
	Line_min = np.zeros(0)
	Line_max = np.zeros(0)
	Line_step = (waves[1]-waves[0])*1e3 # in mA
	
	for i in range(range_wave.shape[0]):
		Line_min  = np.append(Line_min,waves[np.argmin(np.abs(waves-range_wave[i,0]))])
		Line_max  = np.append(Line_max,waves[np.argmin(np.abs(waves-range_wave[i,1]))])
	
	# Define wavelength grid to be saved
	with open(filename, 'w') as f:
		for i in range(len(atoms)):
			ind = np.where(line['Line'] == int(atoms[i].split(',')[0]))[0][0] # Which index in line file corresponds to the atom
			llambdas = (np.array([Line_min[i],Line_max[i]]) - line['wavelength'][ind])*1e3 # in mA ; Determine relative wavelengths
			f.write(f"{atoms[i]}: {'%6.4f' % llambdas[0]},     {'%2.6f' % Line_step},     {'%6.4f' % llambdas[-1]}\n")
	

def write_grid_mc(conf, filename = 'Grid.grid'):
	"""
	Writes the Grid file with data from the config file

	Parameter
	---------
	config : dict
		Dictionary containing all the information from the config file
	filename : string, optional
		String containing the name of the Grid file. Default: Grid.grid
	"""

	# Load data from config
	range_wave = conf['range_wave']
	atoms = conf["atoms"]


	# Define minimum, step and maximum
	Line_min  = np.array(range_wave[:,0]).astype(np.float32)
	Line_step = np.array(range_wave[:,1]).astype(np.float32)
	Line_max  = np.array(range_wave[:,2]).astype(np.float32)
	
	# Define wavelength grid to be saved
	with open(filename, 'w') as f:
		for i in range(len(atoms)):
			f.write(f"{atoms[i]}: {'%6.4f' % Line_min[i]},     {'%2.4f' % Line_step[i]},     {'%6.4f' % Line_max[i]}\n")


def write_control_1c(filename, conf):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameter
	---------
	filename : string
		Save filename of the control file. Typically it is inv.trol
	config : dict
		Dictionary with the information from the config file
	
	"""
	model		= conf['model']		# Base Model
	psf			= conf['psf']		# Use psf by spectral veil
	cycles		= conf['cycles']		# Number of cycles
	weights		= conf['weights']		# Weights in the control file
	nodes_temp	= conf['nodes_temp']	# Nodes in T
	nodes_magn	= conf['nodes_magn']	# Nodes in B
	nodes_vlos	= conf['nodes_vlos']	# Nodes in vlos
	nodes_gamma	= conf['nodes_gamma']	# Nodes in gamma
	nodes_phi		= conf['nodes_phi']		# Nodes in phi
	mu_cos		= conf['mu_cos']		# mu = cos theta
	abundance		= conf['abundance']		# Abundance file
	line			= conf['line']			# Name of the line file
	gas_pressure   = conf['gas_pressure']	# Gas Pressure
	# Write lines
	with open(filename, 'w') as f:
		f.write(f'Number of cycles           (*):{cycles}                  ! (0=synthesis)\n')
		f.write('Observed profiles          (*):' + d.profile_obs + '      ! \n')
		f.write('Stray light file              :                   ! (none=no stray light contam)\n')
		f.write('PSF file                      :' + psf + '        ! (none=no convolution with PSF)\n')
		f.write('Wavelength grid file       (s):' + d.Grid + '! (none=automatic selection)\n')
		f.write('Atomic parameters file        :' + line + '    ! (none=DEFAULT LINES file)\n')
		f.write('Abundances file               :'    + abundance + '         ! (none=DEFAULT ABUNDANCES file)\n')
		f.write('Initial guess model 1      (*):'+ model + '      !\n')
		f.write('Initial guess model 2         :\n')
		f.write('Weight for Stokes I           :'   + weights[0] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes Q           :'   + weights[1] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes U           :'   + weights[2] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes V           :'   + weights[3] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('AUTOMATIC SELECT. OF NODES?   :                   ! (DEFAULT=0=no; 1=yes)\n')
		f.write('Nodes for temperature 1       :'  + nodes_temp + '\n')
		f.write('Nodes for electr. press. 1    :                         \n')
		f.write('Nodes for microturb. 1        :                         \n')
		f.write('Nodes for magnetic field 1    :'+ nodes_magn + '\n')
		f.write('Nodes for LOS velocity 1      :'+ nodes_vlos + '\n')
		f.write('Nodes for gamma 1             :'+ nodes_gamma + '\n')
		f.write('Nodes for phi 1               :'    + nodes_phi + '\n')
		f.write('Invert macroturbulence 1?     :                   ! (0 or blank=no, 1=yes)\n')
		f.write('Nodes for temperature 2       :                   \n')
		f.write('Nodes for electr. press. 2    :                   \n')
		f.write('Nodes for microturb. 2        :                   \n')   
		f.write('Nodes for magnetic field 2    :                   \n')
		f.write('Nodes for LOS velocity 2      :                   \n')
		f.write('Nodes for gamma 2             :                   \n')
		f.write('Nodes for phi 2               :                   \n')
		f.write('Invert macroturbulence 2?     :                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert filling factor?        :                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert stray light factor?    :0                   ! (0 or blank=no, 1=yes)\n')
		f.write('mu=cos (theta)                :'  + mu_cos +  '              ! (DEFAULT: mu=1)\n')
		f.write('Estimated S/N for I           :200                ! (DEFAULT: 1000) \n')
		f.write('Continuum contrast            :                    ! (DEFAULT: not used)\n')
		f.write('Tolerance for SVD             :1.e-4              ! (DEFAULT value: 1e-4)\n')
		f.write('Initial diagonal element      :                    ! (DEFAULT value: 1.e-3)\n')
		f.write('Splines/Linear Interpolation  :                    ! (0 or blank=splines, 1=linear)\n')
		f.write('Gas pressure at surface 1     :' + gas_pressure + '              ! (0 or blank=Pe boundary cond.)\n')
		f.write('Gas pressure at surface 2     :                    ! (0 or blank=Pe boundary cond.\n')
		f.write('Magnetic pressure term?       :                    ! (0 or blank=no, 1=yes\n')
		f.write("NLTE Departures filename      :                    ! blanck= LTE (Ej.) depart_6494.dat'\n")


	f.close()


######################################################################

def write_control_2c(filename, conf):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameter
	---------
	filename : string
		Save filename of the control file. Typically it is inv.trol
	conf : dict
		Dictionary with the information from the config file
	
	"""
	model1		= d.guess1			# Base Model 1
	model2		= d.guess2			# Base Model 2
	psf			= conf['psf']			# Use psf by spectral veil
	cycles		= conf['cycles']		# Number of cycles
	weights		= conf['weights']		# Weights in the control file
	nodes_temp1	= conf['nodes_temp1']	# Nodes in T
	nodes_magn1	= conf['nodes_magn1']	# Nodes in B
	nodes_vlos1	= conf['nodes_vlos1']	# Nodes in vlos
	nodes_gamma1	= conf['nodes_gamma1']	# Nodes in gamma
	nodes_phi1	= conf['nodes_phi1']		# Nodes in phi
	nodes_temp2	= conf['nodes_temp2']	# Nodes in T
	nodes_magn2	= conf['nodes_magn2']	# Nodes in B
	nodes_vlos2	= conf['nodes_vlos2']	# Nodes in vlos
	nodes_gamma2	= conf['nodes_gamma2']	# Nodes in gamma
	nodes_phi2	= conf['nodes_phi2']		# Nodes in phi

	mu_cos		= conf['mu_cos']		# mu = cos theta
	abundance		= conf['abundance']		# Abundance file
	line			= conf['line']			# Name of the line file
	gas_pressure   = conf['gas_pressure']	# Gas Pressure
	# Write lines
	with open(filename, 'w') as f:
		f.write(f'Number of cycles           (*):{cycles}                  ! (0=synthesis)\n')
		f.write('Observed profiles          (*):' + d.profile_obs + '      ! \n')
		f.write('Stray light file              :                   ! (none=no stray light contam)\n')
		f.write('PSF file                      :' + psf + '        ! (none=no convolution with PSF)\n')
		f.write('Wavelength grid file       (s):' + d.Grid + '! (none=automatic selection)\n')
		f.write('Atomic parameters file        :' + line + '    ! (none=DEFAULT LINES file)\n')
		f.write('Abundances file               :' + abundance + '         ! (none=DEFAULT ABUNDANCES file)\n')
		f.write('Initial guess model 1      (*):' + model1 + ' \n')
		f.write('Initial guess model 2         :' + model2 + ' \n')
		f.write('Weight for Stokes I           :' + weights[0] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes Q           :' + weights[1] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes U           :' + weights[2] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes V           :' + weights[3] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('AUTOMATIC SELECT. OF NODES?   :                   ! (DEFAULT=0=no; 1=yes)\n')
		f.write('Nodes for temperature 1       :' + nodes_temp1 + '\n')
		f.write('Nodes for electr. press. 1    :                         \n')
		f.write('Nodes for microturb. 1        :                         \n')
		f.write('Nodes for magnetic field 1    :' + nodes_magn1 + '\n')
		f.write('Nodes for LOS velocity 1      :' + nodes_vlos1 + '\n')
		f.write('Nodes for gamma 1             :' + nodes_gamma1 + '\n')
		f.write('Nodes for phi 1               :' + nodes_phi1 + '\n')
		f.write('Invert macroturbulence 1?     :                   ! (0 or blank=no, 1=yes)\n')
		f.write('Nodes for temperature 2       :'  + nodes_temp2 + ' \n')
		f.write('Nodes for electr. press. 2    :                   \n')
		f.write('Nodes for microturb. 2        :                   \n')
		f.write('Nodes for magnetic field 2    :' + nodes_magn2 + '\n')
		f.write('Nodes for LOS velocity 2      :' + nodes_vlos2 + '\n')
		f.write('Nodes for gamma 2             :' + nodes_gamma2 + '\n')
		f.write('Nodes for phi 2               :' + nodes_phi2 + '\n')
		f.write('Invert macroturbulence 2?     :                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert filling factor?        :1                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert stray light factor?    :0                   ! (0 or blank=no, 1=yes)\n')
		f.write('mu=cos (theta)                :'  + mu_cos +  '              ! (DEFAULT: mu=1)\n')
		f.write('Estimated S/N for I           :200                ! (DEFAULT: 1000) \n')
		f.write('Continuum contrast            :                    ! (DEFAULT: not used)\n')
		f.write('Tolerance for SVD             :1.e-4              ! (DEFAULT value: 1e-4)\n')
		f.write('Initial diagonal element      :                    ! (DEFAULT value: 1.e-3)\n')
		f.write('Splines/Linear Interpolation  :                    ! (0 or blank=splines, 1=linear)\n')
		f.write('Gas pressure at surface 1     :' + gas_pressure + '              ! (0 or blank=Pe boundary cond.)\n')
		f.write('Gas pressure at surface 2     :                    ! (0 or blank=Pe boundary cond.\n')
		f.write('Magnetic pressure term?       :                    ! (0 or blank=no, 1=yes\n')
		f.write("NLTE Departures filename      :                    ! blanck= LTE (Ej.) depart_6494.dat'\n")


	f.close()

######################################################################

def write_control_mc(filename, conf, Type = 'inv'):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameter
	---------
	filename : string
		Save filename of the control file. Typically it is inv.trol
	config : dict
		Dictionary with the information from the config file
	Type : string
		which type of control file is created ('syn' for synthesis, 'inv' for inversion)
	
	"""
	if Type == 'syn':
		model = d.model_syn
	elif Type == 'inv':
		model = d.model_inv
	if Type == 'inv':
		cycles	= conf['cycles']		# Number of cycles
	else:
		cycles 	= 0
	weights		= conf['weights']		# Weights in the control file
	if Type == 'inv':
		nodes_temp	= conf['nodes_temp']	# Nodes in T
		nodes_magn	= conf['nodes_magn']	# Nodes in B
		nodes_vlos	= conf['nodes_vlos']	# Nodes in vlos
		nodes_gamma	= conf['nodes_gamma']	# Nodes in gamma
		nodes_phi		= conf['nodes_phi']		# Nodes in phi
	else:
		nodes_temp	= ''	# Nodes in T
		nodes_magn	= ''	# Nodes in B
		nodes_vlos	= ''	# Nodes in vlos
		nodes_gamma	= ''	# Nodes in gamma
		nodes_phi		= ''		# Nodes in phi

	abundance		= conf['abundance']		# Abundance file
	gas_pressure   = conf['gas_pressure']	# Gas Pressure

	# Write lines
	with open(filename, 'w') as f:
		f.write(f'Number of cycles           (*):{cycles}               ! (0=synthesis)\n')
		f.write('Observed profiles          (*):' + d.profile + '      ! \n')
		f.write('Stray light file              :                   ! (none=no stray light contam)\n')
		f.write('PSF file                      :                   ! (none=no convolution with PSF)\n')
		f.write('Wavelength grid file       (s):' + d.Grid + '     ! (none=automatic selection)\n')
		f.write('Atomic parameters file        :' + conf['line'] + '    ! (none=DEFAULT LINES file)\n')
		f.write('Abundances file               :' + abundance + '         ! (none=DEFAULT ABUNDANCES file)\n')
		f.write('Initial guess model 1      (*):' + model + '      !\n')
		f.write('Initial guess model 2         :\n')
		f.write('Weight for Stokes I           :'   + weights[0] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes Q           :'   + weights[1] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes U           :'   + weights[2] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes V           :'   + weights[3] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('AUTOMATIC SELECT. OF NODES?   :                    ! (DEFAULT=0=no; 1=yes)\n')
		f.write(f'Nodes for temperature 1       :{nodes_temp}\n')
		f.write('Nodes for electr. press. 1    :                         \n')
		f.write('Nodes for microturb. 1        :                         \n')
		f.write(f'Nodes for magnetic field 1    :{nodes_magn}\n')
		f.write(f'Nodes for LOS velocity 1      :{nodes_vlos}\n')
		f.write(f'Nodes for gamma 1             :{nodes_gamma}\n')
		f.write(f'Nodes for phi 1               :{nodes_phi}\n')
		f.write('Invert macroturbulence 1?     :                    ! (0 or blank=no, 1=yes)\n')
		f.write('Nodes for temperature 2       :                    \n')
		f.write('Nodes for electr. press. 2    :                    \n')
		f.write('Nodes for microturb. 2        :                    \n')   
		f.write('Nodes for magnetic field 2    :                    \n')
		f.write('Nodes for LOS velocity 2      :                    \n')
		f.write('Nodes for gamma 2             :                    \n')
		f.write('Nodes for phi 2               :                    \n')
		f.write('Invert macroturbulence 2?     :                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert filling factor?        :                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert stray light factor?    :                    ! (0 or blank=no, 1=yes)\n')
		f.write('mu=cos (theta)                :                    ! (DEFAULT: mu=1)\n')
		f.write('Estimated S/N for I           :200                 ! (DEFAULT: 1000) \n')
		f.write('Continuum contrast            :                    ! (DEFAULT: not used)\n')
		f.write('Tolerance for SVD             :1.e-6               ! (DEFAULT value: 1e-4)\n')
		f.write('Initial diagonal element      :                    ! (DEFAULT value: 1.e-3)\n')
		f.write('Splines/Linear Interpolation  :                    ! (0 or blank=splines, 1=linear)\n')
		f.write('Gas pressure at surface 1     :' + gas_pressure + '          ! (0 or blank=Pe boundary cond.)\n')
		f.write('Gas pressure at surface 2     :                    ! (0 or blank=Pe boundary cond.\n')
		f.write('Magnetic pressure term?       :                    ! (0 or blank=no, 1=yes\n')
		f.write("NLTE Departures filename      :                    ! blanck= LTE (Ej.) depart_6494.dat'\n")


	f.close()



######################################################################
def write_model(filename, Header, log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z = None, Pg = None, rho = None):
	"""
	Write a model with the given data in a specific format. Note that negative values
	have one white space less

	Parameter
	---------
	filename : string
		Name of the saved file
	Header : string
		Header of the model
	log_tau : numpy.array
		Log tau
	T : numpy.array
		Temperature in K
	Pe : numpy.array
		Electron pressure in dyn/cm^2
	v_micro : numpy.array
		Microturbulence velocity in cm/s
	B : numpy.array
		Magnetic field strength in Gauss
	vlos : numpy.array
		Line-of-sight velocity in cm/s
	inc : numpy.array
		Inclination in deg
	azimuth : numpy.array
		Azimuth angle in deg
	z : numpy.array, optional
		Height in km
	Pg : numpy.array, optional
		Gas pressure in dyn/cm^2
	rho : numpy.array, optional
		Density in g/cm^3

	Return
	------
	None
	"""
	if z is not None:
		f = open(filename, 'w')
		f.write(f"{Header}\n")
		for n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11 in zip(log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho):
			f.write(f" {n1:>7.4f} {n2:>7.1f} {n3:>12.5E} {n4:>10.3E} {n5:>11.4E} {n6:>11.4E} {n7:>11.4E} {n8:>11.4E} {n9:>11.4E} {n10:>11.4E} {n11:>11.4E}\n")

	else:
		f = open(filename, 'w')
		f.write(f"{Header}\n")
		for n1,n2,n3,n4,n5,n6,n7,n8 in zip(log_tau, T, Pe, v_micro, B, vlos, inc, azimuth):
			f.write(f" {n1:>7.4f} {n2:>7.1f} {n3:>12.5E} {n4:>10.3E} {n5:>11.4E} {n6:>11.4E} {n7:>11.4E} {n8:>11.4E}\n")

######################################################################

def option(text1, text2):
	"""
	Print an option in a help page

	Parameter
	---------
	text1 : str
		First text
	text2 : str
		Second text
	
	Return
	------
	None
	"""
	print(f"{text1}")
	print(f"\t{text2}")
	

###############################################################################################################################################################

def read_chi2s(conf, tasks):
	"""
	Reads all the chi2 from the inversion
	
	Parameter
	---------
	config : dict
		Config parameters
	tasks : dict
		Dictionary with the used folders
	
	Return
	------
	chi2 : numpy array
		Numpy array with all chi2 values
	
	"""
	if conf['mode'] != "MC":
		print("[ERROR] This function 'read_chi2s' is only defined for the mode MC")
		return np.empty()
	
	filename = d.inv_trol_file[:d.inv_trol_file.rfind('.')] + ".chi"
	path = conf['path']
	num = conf['num']
	

	chi2 = np.zeros(shape=(num))
	for i in range(num):
		chi2[i] = read_chi2(f"{os.path.join(path,tasks['folders'][i])}/{filename}")

	return chi2

######################################################################

def write_profile(filename, profiles, pos):
	"""
	Write a profile for a specific model number to a file

	Parameter
	---------
	filename : string
		Name of the saved file
	profiles : list
		List containing all the profiles
	atoms : list
		List containing the number of the line from the Line file
	pos : int
		Position which model is saved

	Return
	------
	None
	"""

	f = open(filename, 'w')
	#for a in range(len(atoms)):
	for m in range(profiles.shape[2]):
			f.write(f" {int(profiles[pos,0,m]):>2} {profiles[pos,1,m]:>10.4f} {profiles[pos,2,m]:>14.7E} {profiles[pos,3,m]:>14.7E} {profiles[pos,4,m]:>14.7E} {profiles[pos,5,m]:>14.7E}\n")
	f.close()

######################################################################

def option(text1, text2):
	"""
	Print an option in a help page

	Parameter
	---------
	text1 : str
		First text
	text2 : str
		Second text
	
	Return
	------
	None
	"""
	print(f"{text1}")
	print(f"\t{text2}")


######################################################################################################################### 2 components



def read_chi2(filename, task = ''):
	"""
	Reads the last chi value in a inv.chi file
	
	Parameter
	---------
	filename : string
		Path of the chi file
	task : string, optional
		Prints out in which folder the chi2 file does not exist. Default: ''

	
	Return
	------
	chi2 : float
		Best chi2 value of the fit
	
	
	"""
	if not exists(filename):
		print("[read_chi2] " + filename + " does not exist in " + task + ".")
		sys.exit(1)

	# Load data
	data = np.genfromtxt(filename)
	return data[-1][1]
