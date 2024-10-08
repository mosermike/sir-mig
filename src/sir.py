"""
Miscellaneous Functions used throughout the code
"""
import numpy as np 
import sys
import os
from os.path import exists
from typing import Tuple

def create_task_folder_list(arg):
	"""
	Creates a list which folders should be created and executed. This is done so
	that the inversion itself can be executed linearly to make use of all cores.

	Parameters
	----------
	arg : numpy array or int
		1x4 array containing the limits in x and y of the data or number of 1D models

	Returns
	-------
	out : dict
		Dictionary with all the names of the task folders, x and y position
	"""
	import numpy as np
	import definitions as d
	# Create arrays
	tasks = []
	xs = []
	ys = []
	
	# Task folder list for mode MC
	if isinstance(arg,int):
		y = 0
		# Determine task folder names
		for x in range(arg):
			x_str, y_str = x_y_add_zeros(x,0)

			tasks.append(d.task_start + x_str + "_" + y_str)
			xs.append(x)
			ys.append(y)
		Dict = {
				'folders' : tasks,
				'x' : np.array(xs),
				'y' : np.array(ys),
		
		}

	# Task List for mode 1C and 2C
	else:
		Map = arg
		if Map[1] < Map[0]:
			raise Exception(f"xmax in the map argument is smaller than xmin")
		if Map[3] < Map[2]:
			raise Exception(f"ymax in the map argument is smaller than ymin")
		# Determine task folder names
		for x in range(Map[0], Map[1]+1):
			for y in range(Map[2], Map[3]+1):
				x_str, y_str = x_y_add_zeros(x, y)

				tasks.append(d.task_start + x_str + "_" + y_str)
				xs.append(x)
				ys.append(y)
		Dict = {
				'folders': tasks,
				'x': np.array(xs),
				'y': np.array(ys),
		
		}

	return Dict

def determine_line_core(linefile : str, num : int) -> float:
	"""
	Determines the spectral line core value from a number in the line file

	Parameters
	----------
	linefile : str
		Path to the linefile
	num : int
		Number of the line core

	Returns
	-------
	determine_line_core : float
		Spectral Core Number of the desired number
	"""

	line = read_line(linefile)
	for i in range(len(line["Line"])):
		if line["Line"][i] == num:
			return line["wavelength"][i]
	print(f"[determine_line_core] The number {num} does not exist in the provided lines file {linefile}")
	return 0

def initial(mode : str):
	"""
	Initial print outs and preparation

	Parameters
	----------
	mode : str
		Mode which is used
	
	Returns
	-------
	None

	

	"""
	print()
	print("╭───────────────────────────────────────────────────╮")
	print("│ SIR - MIG                                         │")
	print("│ Version 1.0                                       │")
	print("│ Multiple Initial Guesses                          │")
	print(f"│ Mode: {mode}                                          │")
	print("│ Author: Mike Moser                                │")
	print("╰───────────────────────────────────────────────────╯")	
	print()

def mpl_library():
	"""
	Adjust the matplotlib settings to the defined library in the definitions file or to the default one

	

	"""
	import matplotlib.pyplot as plt
	import definitions as d
	
	dirname = os.path.dirname(os.path.realpath(__file__))
	plt.rcParams["savefig.format"] = "pdf"
	if d.plt_lib != "":
		plt.style.use(d.plt_lib)
	else:
		if os.path.exists(dirname + '/mml.mplstyle'):
			plt.style.use(dirname + '/mml.mplstyle')
			# if dvipng is not installed, dont use latex
			import shutil
			if shutil.which('dvipng') is None:
				plt.rcParams["text.usetex"] = "False"
				plt.rcParams["font.family"] = 'sans-serif'
				plt.rcParams["mathtext.fontset"] = 'dejavuserif'
		elif "mml" in plt.style.available:
			plt.style.use('mml')
			# if dvipng is not installed, dont use latex
			import shutil
			if shutil.which('dvipng') is None:
				plt.rcParams["text.usetex"] = "False"
				plt.rcParams["font.family"] = 'sans-serif'
				plt.rcParams["mathtext.fontset"] = 'dejavuserif'
	return

def option(text1 : str, text2 : str) -> None:
	"""
	Print an option in a help page

	Parameters
	----------
	text1 : str
		First text
	text2 : str
		Second text
	
	Returns
	-------
	None

	

	"""
	print(f"{text1}")
	print(f"\t{text2}")

def read_chi2(filename, task : str= '') -> float:
	"""
	Reads the last chi value in a inv.chi file
	
	Parameters
	----------
	filename : string
		Path of the chi file
	task : string, optional
		Prints out in which folder the chi2 file does not exist, by default ''

	
	Returns
	-------
 	float
		Best chi2 value of the fit
	
	
	
	"""
	
	if not exists(filename):
		print("[read_chi2] " + filename + " does not exist in " + task + ".")
		return -1

	# Load data
	data = np.genfromtxt(filename)
	return data[-1][1]


def read_config(filename : str, check : bool = False) -> dict:
	"""
	Reads a config file for the inversion
	
	Parameters
	----------
	filename : string
		Path of the control file
	check : bool, optional
		Check if file exists, by default 'False'

	Returns
	-------
	dict
		Dict containing all the information from the config file
	
	Raises
	------
	FileExistsError
		if config file does not exist
	FileNotFoundError
		if path does not exist
	Exception
		if range_wave and atoms do not consist of the same number of lines
	FileNotFoundError
		if the synthesis input does not exist in mode SY
	
	"""
	if not exists(filename):
		raise FileExistsError("[read_config] " + filename + " does not exist.")

	# Load data
	data = np.genfromtxt(filename, delimiter=':', comments='#', dtype=str)

	# Remove spaces in first column
	for i in range(len(data)):
		# Check if the value is empty
		if len(data[i][0].replace(' ','')) == 0:
			data[i][0] = ''
			continue
		while data[i][0][0] == ' ':
			data[i][0] = data[i][0][1:]
		while data[i][0][-1] == ' ':
			data[i][0] = data[i][0][:-1]

	# Remove spaces in second column
	for i in range(len(data)):
		# Check if the value is empty
		if len(data[i][1].replace(' ','')) == 0:
			data[i][1] = ''
			continue
		while data[i][1][0] == ' ':
			data[i][1] = data[i][1][1:]
		while data[i][1][-1] == ' ':
			data[i][1] = data[i][1][:-1]

		data[i][0] = data[i][0].replace('  ','')   # replace double spaces

	# Create dictionary
	Dict = {
		"filename" : filename,
	}

	for i in data:
		if i[0] != '':
			Dict[i[0]] = i[1]

	# replace path if it is ./ but there is a path in the filename
	if Dict["path"][:2] == "./" and "/" in filename:
		Dict["path"] = os.path.join(filename[:filename.rfind("/")], Dict["path"][2:])

	if check:
		if(Dict["mode"] != "SY"):
			if Dict["inv_out"] == "":
				print("[read_config] No value in inv_out")
			if Dict["cycles"] == '':
				print("[read_config] No value in cycles")
	if Dict["line"] == "":
		print("[read_config] No value in line")
	if Dict["atoms"] == "":
		print("[read_config] No value in atoms")
	if check:
		if Dict["abundance"] == '':
			print("[read_config] No value in abundance")

	# Transform the information into lists or in different types than string
	if Dict['mode'] == "1C" or Dict['mode'] == "2C":
		Dict['map'] = np.array([int(i) for i in Dict["map"].split(',')], dtype=int)
		
		if "quiet_sun" in Dict:
			Dict["quiet_sun"] = np.array([int(i) for i in Dict["quiet_sun"].split(',')])
		
	if Dict['mode'] == "MC":
		Dict['num'] = int(Dict['num'])

	temp = [i. split(',') for i in Dict["range_wave"].split(';')]
	Dict["range_wave"] = np.zeros((len(temp),3), dtype=np.float64)
	for i in range(0,len(temp)):
		Dict["range_wave"][i,0] = np.float64(temp[i][0])
		Dict["range_wave"][i,1] = np.float64(temp[i][1])
		Dict["range_wave"][i,2] = int(temp[i][2].replace(".0",""))

	Dict['atoms'] = Dict["atoms"].split(';')	# Atoms

	if(Dict["mode"] != "SY"):
		Dict["random_guess"] = int(Dict["random_guess"])

		Dict["random_pars"] = Dict["random_pars"].replace(" ",'').split(",")
	
		Dict["cycles"] = int(Dict["cycles"])

	Dict["weights"] = Dict["weights"].split(',')

	# TODO check the checks
	if check:
		if len(Dict['atoms']) != len(Dict['range_wave']):
			raise Exception("[read_config] The number of lines in 'atoms' do not fit the given ranges in 'range_wave'! Abort...")
		if not exists(Dict['path']):
			raise FileNotFoundError(f"[read_config] {Dict['path']} does not exist.")
		if (Dict["mode"] == "1C" or Dict["mode"] == "MC") and not exists(os.path.join(Dict['path'],Dict['model'])):
			print(f"[read_config] {Dict['model']} does not exist.")
		if Dict["mode"] == "2C" and not exists(os.path.join(Dict['path'],Dict['model1'])):
			print(f"[read_config] {Dict['model1']} does not exist.")
		if Dict["mode"] == "2C" and not exists(os.path.join(Dict['path'],Dict['model2'])):
			print(f"[read_config] {Dict['model2']} does not exist.")
		if Dict["mode"] == "SY" and not exists(os.path.join(Dict['path'],Dict['syn_in'])):
			raise FileNotFoundError(f"File '{os.path.join(Dict['path'],Dict['syn_in'])}' does not exist")
		if Dict['mode'] == "1C" or Dict['mode'] == "2C":
			if not exists(os.path.join(Dict['path'],Dict['cube'])) and Dict['preprocess'] == "0":
				raise FileNotFoundError(f"File '{os.path.join(Dict['path'],Dict['syn_in'])}' does not exist")
			if Dict['preprocess'] == "1" and not exists(os.path.join(Dict['path'],Dict['fts_file'])):
				print(f"[read_config] {Dict['fts_file']} does not exist.")
		if Dict['mode'] != "SY":
			if (Dict['chi2'] != "" and Dict['chi2'] != "0" and Dict['chi2'] != "1"):
				print(f"[read_config] Unknown option '{Dict['chi2']}'. chi2 is not computed.")

		if Dict["mode"] == "1C" or Dict["mode"] == "MC":
			if len(Dict["nodes_temp"].split(",")) > Dict["cycles"] or len(Dict["nodes_magn"].split(",")) > Dict["cycles"] or len(Dict["nodes_vlos"].split(",")) > Dict["cycles"] or len(Dict["nodes_gamma"].split(",")) > Dict["cycles"] or len(Dict["nodes_phi"].split(",")) > Dict["cycles"]:
				print("[read_config] Warning: More nodes specified than number of cycles.")
		if Dict["mode"] == "2C":
			if len(Dict["nodes_temp1"].split(",")) > Dict["cycles"] or len(Dict["nodes_magn1"].split(",")) > Dict["cycles"] or len(Dict["nodes_vlos1"].split(",")) > Dict["cycles"] or len(Dict["nodes_gamma1"].split(",")) > Dict["cycles"] or len(Dict["nodes_phi1"].split(",")) > Dict["cycles"]:
				print("[read_config] Warning: More nodes specified for model 1 than number of cycles.")
			if len(Dict["nodes_temp2"].split(",")) > Dict["cycles"] or len(Dict["nodes_magn2"].split(",")) > Dict["cycles"] or len(Dict["nodes_vlos2"].split(",")) > Dict["cycles"] or len(Dict["nodes_gamma2"].split(",")) > Dict["cycles"] or len(Dict["nodes_phi2"].split(",")) > Dict["cycles"]:
				print("[read_config] Warning: More nodes specified for model 2 than number of cycles.")
	

	return Dict


def read_control(filename : str) -> dict:
	"""
	Reads a control file in the scheme SIR expects it.
	
	Parameters
	----------
	filename : str
		Path of the control file
	
	Returns
	-------
	dict
		Dict containing all the information from the control file
	
	Raises
	------
	FileExistsError
		if file does not exist

	

	"""
	if not exists(filename):
		raise FileExistsError("[read_control] " + filename + " does not exist.")

	# Load data
	data = np.genfromtxt(filename, delimiter=':', comments='!', dtype=str)
	for i in range(len(data)):
		data[i][0] = data[i][0].replace('(*)','') # replace the (*)
		data[i][0] = data[i][0].replace('  ','') # replace double spaces
	
	Dict = dict()
	for i in data:
		Dict[i[0]] = i[1]
	return Dict

def read_info(filename : str) -> dict:
	"""
	Reads the info file created in the preprocess

	Parameters
	----------
	filename : str
		File name to be loaded

	Returns
	-------
	dict
		Dictionary with the informations from the info file

	Raises
	------
	FileExistsError
		if file does not exist

	
	"""
	if not exists(filename):
		raise FileExistsError("[read_info] " + filename + " does not exist.")
	infos = dict(np.genfromtxt(filename, dtype='str', delimiter="="), dtype=str)

	return infos
	
def read_model(filename : str):
	"""
	Reads a model file and returns all parameters
	
	Parameters
	----------
	filename : str
		String containing the path of the file

	Returns
	-------
	numpy.array
		Log tau
	numpy.array
		Temperature in K
	numpy.array
		Electron pressure in dyn/cm^2
	numpy.array
		Microturbulence velocity in cm/s
	numpy.array
		Magnetic field strength in Gauss
	numpy.array
		Line-of-sight velocity in cm/s
	numpy.array
		Inclination in deg
	numpy.array
		Azimuth angle in deg
	numpy.array, optional
		Height in km
	numpy.array, optional
		Gas pressure in dyn/cm^2
	numpy.array, optional
		Density in g/cm^3

	Raises
	------
	FileExistsError
		if file does not exist

	
	"""	
	if not exists(filename):
		raise FileExistsError("[read_model] " + filename + " does not exist.")
	
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


def read_profile(filename : str, num : int = 0):
	"""
	Reads the first LINE data from a profile computed by SIR
	
	Parameters
	----------
	filename : str
		String containing the path of the file
	num : int, optional
		Number of the line which is loaded, by default 0 (use first one from line)

	Returns
	-------
	numpy.array
		Wavelengths in A
	numpy.array
		Stokes I
	numpy.array
		Stokes Q
	numpy.array
		Stokes U
	numpy.array 
		Stokes V
	
	Raises
	------
	FileExistsError
		if file does not exist

	
	"""
	if not exists(filename):
		raise FileExistsError("[read_profile] " + filename + " does not exist.")
	
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


def read_grid(filename : str) -> dict:
	"""
	Reads the grid file
	
	Parameters
	----------
	filename : str
		File to be read
	
	Returns
	-------
	dict
		Dict. with 'Line', 'min', 'step' and 'max' in it

	Raises
	------
	FileExistsError
		if file does not exist

	
	"""
	if not exists(filename):
		raise FileExistsError("[read_grid] " + filename + " does not exist.")
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


def read_line(filename):
	"""
	Reads the line file

	Parameters
	----------
	filename : string
		File to be read
	
	Returns
	-------
	dict
		Dict. with 'Line', 'Ion', 'wavelength', 'factor', 'Exc_Pot', log_gf',
		'Transition', 'alpha' and 'sigma' in it
	
	Raises
	------
	FileExistsError
		if file does not exist

	
	"""
	if not exists(filename):
		raise FileExistsError("[read_line] " + filename + " does not exist.")
	
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


def list_to_string(temp, let = ','):
	"""
	Convert a list to a string

	Parameters
	----------
	temp : list
		List with the information which is converted into a string
	let : str
		Letter which is added as a separation

	Returns
	-------
	str
		Information from the list

	

	"""
	if(isinstance(temp,str)):
		return temp
	temp1 = ''
	for i in range(len(temp)):
		temp1 += str(temp[i])
		if i + 1 < len(temp):
				temp1 += let
	return temp1


def _write_config_1c(File : str, conf : dict, verbose : bool = True):
	"""
	Writes a config file with the information provided as a dictionary for the mode 1C

	Parameters
	----------
	File : str
		Save path
	conf : dict
		Dictionary with all the informations
	verbose : bool,optional
		Verbose output that config is written and possibility to abort, by default True

	
	"""
	if verbose:
		print("[write_config] Manually added comments will be overwritten? 1s to abort left ...")
		import time
		time.sleep(2)
	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	if(isinstance(temp,str)):
		range_wave = temp
	else:
		for i in range(len(temp)):
			for j in range(len(temp[i])):
				# Last entry is an integer
				if j == len(temp[i])-1:
					range_wave += str(int(temp[i][j]))
				else:
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
		f.write("# SIR MIG config file\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")

		f.write(f"# \n")
		f.write(f"# Stuff from the data\n")
		f.write(f"# \n")
		f.write(f"cube : {conf['cube']} # Data cube used for the inversion (bin)\n")
		f.write(f"map  : {Map} # Pixels to be considered as a list\n")

		f.write("#\n")
		f.write("# Data Preprocessing\n")
		f.write("#\n")
		f.write(f"preprocess : {conf['preprocess']} # Preprocess data (1 = True, 0 = False)\n")
		f.write(f"instrument : {conf['instrument']} # Instrument used (GRIS, Hinode or empty)\n")
		f.write(f"ending     : {conf['ending']} # Ending of GRIS file used for merging\n")
		f.write(f"quiet_sun  : {quiet_sun} # Quiet sun region for normalization as a list (0 => already normalised)\n")
		f.write(f"fts_file   : {conf['fts_file']} # Absolute path to fts file, blank = do not correct spectral veil\n")
		f.write(f"shift_wave : {conf['shift_wave']} # Shift the wavelength grid when waves file is created in mA\n")
		f.write(f"save_cube  : {conf['save_cube']} # Save preprocessed data (merged and normalised) (1=True,0=False)\n")

		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"model      : {conf['model']} # Base Model for guess\n")
		f.write(f"range_wave : {range_wave} # Range for the grid file as (Start wavelength in abs. wavelength, Step in mA, Number of wavelenghts) for each line in the grid file.\n")
		f.write(f"inv_out    : {conf['inv_out']} # Prefix of output of the inversion files\n")
		f.write(f"chi2       : {conf['chi2']} # Compute chi2 and save it under this name\n")
		f.write(f"line       : {conf['line']} # Line file\n")
		f.write(f"atoms      : {atoms} # Atoms used, ; defines a new line\n")
		f.write(f"guess      : {conf['guess']} # Use a bin file as initial guesses, blank use base model\n")
		f.write(f"psf        : {conf['psf']} # Spectral PSF .dat file, 'gauss 1.0' or blank=not used\n")

		f.write(f"# \n")
		f.write(f"# Control file\n")
		f.write(f"# \n")
		f.write(f"cycles       : {conf['cycles']} # Number of cycles\n")
		f.write(f"weights      : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp   : {conf['nodes_temp']} # Nodes in T\n")
		f.write(f"nodes_magn   : {conf['nodes_magn']} # Nodes in B\n")
		f.write(f"nodes_vlos   : {conf['nodes_vlos']} # Nodes in vlos\n")
		f.write(f"nodes_gamma  : {conf['nodes_gamma']} # Nodes in gamma\n")
		f.write(f"nodes_phi    : {conf['nodes_phi']} # Nodes in phi\n")
		f.write(f"vmacro       : {conf['vmacro']} # Macroturbulence velocity\n")
		f.write(f"mu_cos       : {conf['mu_cos']} # mu = cos theta\n")
		f.write(f"abundance    : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition\n")

		f.write(f"# \n")
		f.write(f"# Radomisation Settings\n")
		f.write(f"# \n")
		f.write(f"random_guess : {conf['random_guess']} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars  : {random_pars} # Randomise these parameters in the file(s) below\n")
		f.write(f"lim_B        : {conf['lim_B']} # Limits for the randomisation of B in G\n")
		f.write(f"lim_vlos     : {conf['lim_vlos']} # Limits for the randomisation of vlos in cm/s\n")
		f.write(f"lim_gamma    : {conf['lim_gamma']} # Limits for the randomisation of the inclination in deg\n")
		f.write(f"lim_phi      : {conf['lim_phi']} # Limits for the randomisation of the azimuth in deg")

def _write_config_2c(File : str, conf : dict, verbose : bool = True):
	"""
	Writes a config file with the information provided as a dictionary

	Parameters
	----------
	File : str
		Save path
	conf : dict
		Dictionary with all the informations
	verbose : bool,optional
		Verbose output that config is written and possibility to abort, by default True

	

	"""
	if verbose:
		print("[write_config] Manually added comments will be overwritten? 1s to abort left ...")
		import time
		time.sleep(2)
	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	if(isinstance(temp,str)):
		range_wave = temp
	else:
		for i in range(len(temp)):
			for j in range(len(temp[i])):
				# Last entry is an integer
				if j == len(temp[i])-1:
					range_wave += str(int(temp[i][j]))
				else:
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
		f.write("# SIR MIG config file\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")

		f.write(f"# \n")
		f.write(f"# Stuff from the data\n")
		f.write(f"# \n")
		f.write(f"cube : {conf['cube']} # Data cube used for the inversion (bin)\n")
		f.write(f"map  : {Map} # Pixels to be considered as a list (0 means all pixels)\n")

		f.write("#\n")
		f.write("# Data Preprocessing\n")
		f.write("#\n")
		f.write(f"preprocess : {conf['preprocess']} # Preprocess data (1 = True, 0 = False)\n")
		f.write(f"instrument : {conf['instrument']} # Instrument used (GRIS, Hinode or empty)\n")
		f.write(f"ending     : {conf['ending']} # Ending of GRIS file used for merging\n")
		f.write(f"quiet_sun  : {quiet_sun} # Quiet sun region for normalization as a list (0 => already normalised)\n")
		f.write(f"fts_file   : {conf['fts_file']} # Absolute path to fts file, blank = do not correct spectral veil\n")
		f.write(f"shift_wave : {conf['shift_wave']} # Shift the wavelength grid when waves file is created in mA\n")
		f.write(f"save_cube  : {conf['save_cube']} # Save preprocessed data (merged and normalised) (1=True,0=False)\n")

		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"model1     : {conf['model1']} # Base Model 1 for guess\n")
		f.write(f"model2     : {conf['model2']} # Base Model 2 for guess\n")
		f.write(f"range_wave : {range_wave} # Range for the grid file as (Start wavelength in abs. wavelength, Step in mA, Number of wavelenghts) for each line in the grid file.\n")
		f.write(f"fill       : {conf['fill']} # Filling factors for both models separated by a ',' (if random_guess > 0)\n")
		f.write(f"inv_out    : {conf['inv_out']} # Prefix of output of the inversion files\n")
		f.write(f"chi2       : {conf['chi2']} # Compute chi2 and save it under this name\n")
		f.write(f"line       : {conf['line']} # Line file\n")
		f.write(f"atoms      : {atoms} # Atoms used, ; defines a new line\n")
		f.write(f"guess1     : {conf['guess1']} # Use a bin file as initial guesses, blank use base model 1\n")
		f.write(f"guess2     : {conf['guess2']} # Use a bin file as initial guesses, blank use base model 2\n")
		f.write(f"psf        : {conf['psf']} # Spectral PSF .dat file, 'gauss 1.0' or blank=not used\n")

		f.write(f"# \n")
		f.write(f"# Control file\n")
		f.write(f"# \n")
		f.write(f"cycles       : {conf['cycles']} # Number of cycles\n")
		f.write(f"weights      : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp1  : {conf['nodes_temp1']} # Nodes 1 in T\n")
		f.write(f"nodes_magn1  : {conf['nodes_magn1']} # Nodes 1 in B\n")
		f.write(f"nodes_vlos1  : {conf['nodes_vlos1']} # Nodes 1 in vlos\n")
		f.write(f"nodes_gamma1 : {conf['nodes_gamma1']} # Nodes 1 in gamma\n")
		f.write(f"nodes_phi1   : {conf['nodes_phi1']} # Nodes in 1 phi\n")
		f.write(f"nodes_temp2  : {conf['nodes_temp2']} # Nodes 2 in T\n")
		f.write(f"nodes_magn2  : {conf['nodes_magn2']} # Nodes 2 in B\n")
		f.write(f"nodes_vlos2  : {conf['nodes_vlos2']} # Nodes 2 in vlos\n")
		f.write(f"nodes_gamma2 : {conf['nodes_gamma2']} # Nodes 2 in gamma\n")
		f.write(f"nodes_phi2   : {conf['nodes_phi2']} # Nodes 2 in phi\n")
		f.write(f"invert_fill  : {conf['invert_fill']} # Invert filling factor (1 or 0)")
		f.write(f"mu_cos       : {conf['mu_cos']} # mu = cos theta\n")
		f.write(f"abundance    : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition (two entries separated with ' ')\n")

		f.write(f"# \n")
		f.write(f"# Radomisation Settings\n")
		f.write(f"# \n")
		f.write(f"random_guess : {conf['random_guess']} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars  : {random_pars} # Randomise these parameters in the file(s) below\n")
		f.write(f"lim_B1       : {conf['lim_B1']} # Limits 1 for the randomisation of B in G\n")
		f.write(f"lim_vlos1    : {conf['lim_vlos1']} # Limits 1 for the randomisation of vlos in cm/s\n")
		f.write(f"lim_gamma1   : {conf['lim_gamma1']} # Limits 1 for the randomisation of the inclination in deg\n")
		f.write(f"lim_phi1     : {conf['lim_phi1']} # Limits 1 for the randomisation of the azimuth in deg\n")
		f.write(f"lim_B2       : {conf['lim_B2']} # Limits 2 for the randomisation of B in G\n")
		f.write(f"lim_vlos2    : {conf['lim_vlos2']} # Limits 2 for the randomisation of vlos in cm/s\n")
		f.write(f"lim_gamma2   : {conf['lim_gamma2']} # Limits 2 for the randomisation of the inclination in deg\n")
		f.write(f"lim_phi2     : {conf['lim_phi2']} # Limits 2 for the randomisation of the azimuth in deg")

def _write_config_mc(File : str, conf : dict, verbose : bool=True):
	"""
	Writes a config file with the information provided as a dictionary

	Parameters
	----------
	File : str
		Save path
	conf : dict
		Dictionary with all the informations
	verbose : bool,optional
		Verbose output that config is written and possibility to abort

	

	"""
	if verbose:
		print("[write_config] Note that manually added comments will be overwritten! 1s left to abort ...")
		import time
		time.sleep(2)
	
	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	if(isinstance(temp,str)):
		range_wave = temp
	else:
		for i in range(len(temp)):
			for j in range(len(temp[i])):
				# Last entry is an integer
				if j == len(temp[i])-1:
					range_wave += str(int(temp[i][j]))
				else:
					range_wave += str(temp[i][j])
				if j +1 < len(temp[i]):
					range_wave += ','
			if i < len(temp)-1:
				range_wave += ';'

	atoms = list_to_string(conf["atoms"], ";")
	weights = list_to_string(conf["weights"])
	random_pars = list_to_string(conf["random_pars"])
	
	with open(File, 'w') as f:
		f.write(f"# SIR MIG config file\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")

		f.write(f"# \n")
		f.write(f"# General Stuff\n")
		f.write(f"# \n")
		f.write(f"num        : {conf['num']} # Number of Models\n")
		f.write(f"model      : {conf['model']} # Base Model for guess\n")
		f.write(f"atoms      : {atoms} # Atoms to be used in Grid file\n")	
		f.write(f"range_wave : {range_wave} # Ranges of wavelengths in mA to be considered start1,step1,num1;start2,step2,num2;... First pair belongs to first line in Grid file, etc.\n")
		
		f.write(f"# \n")
		f.write(f"# Data Stuff\n")
		f.write(f"# \n")
		f.write(f"syn_out   : {conf['syn_out']} # Output prefix of the synthesis profiles and models\n")
		f.write(f"noise_out : {conf['noise_out']} # Output prefix of the noise profiles\n")
		f.write(f"inv_out   : {conf['inv_out']} # Prefix of the output of the inversion files\n")
		f.write(f"chi2      : {conf['chi2']} # Compute chi2 and save it under this name\n")

		f.write(f"#\n")
		f.write(f"# Creating Models and Synthesis\n")
		f.write(f"#\n")
		f.write(f"model_nodes   : {conf['model_nodes']} # Create models with 1, 2 or 3 nodes\n")
		f.write(f"model_pars    : {conf['model_pars']} # Randomise these parameters while creating models as a list\n")
		f.write(f"noise_I       : {conf['noise_I']} # Noise in I\n")		
		f.write(f"noise_Q       : {conf['noise_Q']} # Noise in Q\n")
		f.write(f"noise_U       : {conf['noise_U']} # Noise in U\n")
		f.write(f"noise_V       : {conf['noise_V']} # Noise in V\n")
		f.write(f"create_B      : {conf['create_B']} # The limits for the first and last node in B\n")
		f.write(f"create_vlos   : {conf['create_vlos']} # The limits for the first and last node in vlos\n")
		f.write(f"create_gamma  : {conf['create_gamma']} # The limits for the first and last node in gamma\n")
		f.write(f"create_phi    : {conf['create_phi']} # The limits for the first and last node in phi\n")
		f.write(f"create_points : {conf['create_points']} # At this log tau points the models are interpolated with splines (decreasing), 2 or 3 values for 2 or 3 nodes\n")

		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"line         : {conf['line']} # Line file\n")
		f.write(f"guess        : {conf['guess']} # Use a bin file as initial guesses, blank use base model\n")
		f.write(f"cycles       : {conf['cycles']} # Number of cycles\n")
		f.write(f"weights      : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp   : {conf['nodes_temp']} # Nodes in T\n")
		f.write(f"nodes_magn   : {conf['nodes_magn']} # Nodes in B\n")
		f.write(f"nodes_vlos   : {conf['nodes_vlos']} # Nodes in vlos\n")
		f.write(f"nodes_gamma  : {conf['nodes_gamma']} # Nodes in gamma\n")
		f.write(f"nodes_phi    : {conf['nodes_phi']} # Nodes in phi\n")
		f.write(f"vmacro       : {conf['vmacro']} # Macroturbulence velocity\n")
		f.write(f"abundance    : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition\n")

		f.write(f"# \n")
		f.write(f"# Randomisation Settings\n")
		f.write(f"# \n")
		f.write(f"random_guess : {conf['random_guess']} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars  : {random_pars} # Randomise these parameters for the guess\n")
		f.write(f"lim_B        : {conf['lim_B']} # Limits for the randomisation of B in G\n")
		f.write(f"lim_vlos     : {conf['lim_vlos']} # Limits for the randomisation of vlos in cm/s\n")
		f.write(f"lim_gamma    : {conf['lim_gamma']} # Limits for the randomisation of the inclination in deg\n")
		f.write(f"lim_phi      : {conf['lim_phi']} # Limits for the randomisation of the azimuth in deg")

def _write_config_sy(File : str, conf : dict, verbose : bool=True):
	"""
	Writes a config file with the information provided as a dictionary for the mode SY

	Parameters
	----------
	File : str
		Save path
	conf : dict
		Dictionary with all the informations
	verbose : bool,optional
		Verbose output that config is written and possibility to abort

	

	"""
	if verbose:
		print("[write_config] Note that manually added comments will be overwritten! 1s left to abort ...")
		import time
		time.sleep(2)
	
	# Revert, range_wave, map and weights
	range_wave = ''
	temp = conf["range_wave"]
	if(isinstance(temp,str)):
		range_wave = temp
	else:
		for i in range(len(temp)):
			for j in range(len(temp[i])):
				# Last entry is an integer
				if j == len(temp[i])-1:
					range_wave += str(int(temp[i][j]))
				else:
					range_wave += str(temp[i][j])
				if j +1 < len(temp[i]):
					range_wave += ','
			if i < len(temp)-1:
				range_wave += ';'

	atoms = list_to_string(conf["atoms"], ";")
	weights = list_to_string(conf["weights"])

	with open(File, 'w') as f:
		f.write(f"# SIR MIG config file\n")
		f.write(f"mode : {conf['mode']} # Determines which code is executed\n")
		f.write(f"path : {conf['path']} # Path location where all the data is stored and will be saved\n")

		f.write(f"# \n")
		f.write(f"# Data Stuff\n")
		f.write(f"# \n")
		f.write(f"syn_in  : {conf['syn_in']} # Input synthesis models\n")
		f.write(f"syn_out : {conf['syn_out']} # Output of the synthesis profiles\n")
		
		f.write(f"# \n")
		f.write(f"# Synthesis configuration\n")
		f.write(f"# \n")
		f.write(f"atoms        : {atoms} # Atoms to be used in Grid file\n")	
		f.write(f"range_wave   : {range_wave} # Ranges of wavelengths in mA to be considered start1,step1,num1;start2,step2,num2;... First pair belongs to first line in Grid file, etc.\n")
		f.write(f"weights      : {weights} # Weights in the control file\n")
		f.write(f"line         : {conf['line']} # Line file\n")
		f.write(f"vmacro       : {conf['vmacro']} # Macroturbulence velocity\n")
		f.write(f"abundance    : {conf['abundance']} # Abundance file\n")
		f.write(f"gas_pressure : {conf['gas_pressure']} # Gas Pressure Boundary condition")

	return

def write_config(File : str, conf : dict, verbose : bool = True) -> None:
	"""
	Writes a config file with the information provided as a dictionary

	Parameters
	----------
	File : str
		Save path
	conf : dict
		Dictionary with all the informations
	verbose : bool,optional
		Verbose output, by default True.

	Raises
	------
	ValueError
		Mode is not defined
	"""
	if conf["mode"] == "MC":
		_write_config_mc(File, conf, verbose)
	elif conf["mode"] == "1C":
		_write_config_1c(File, conf, verbose)
	elif conf["mode"] == "2C":
		_write_config_2c(File, conf, verbose)
	elif conf["mode"] == "SY":
		_write_config_sy(File, conf, verbose)
	else:
		raise ValueError(f"[write_config] Mode '{conf['mode']}' is not defined and config file cannot be written.")

def write_control(filename : str, conf : dict, Type : str= 'inv'):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameters
	----------
	filename : str
		Save filename of the control file. Typically it is inv.trol
	config : dict
		Dictionary with the information from the config file
	Type : str, optional
		Mode for the sir control file. Options are
		- 'syn': Synthesis
		- 'inv': Inversion (default)
		Only used for mode 'MC'

	Raises
	------
	ValueError
		Mode is not defined

	"""
	if conf['mode'] == 'MC':
		_write_control_mc(filename,conf,Type)
	elif conf['mode'] == '1C':
		_write_control_1c(filename,conf)
	elif conf['mode'] == '2C':
		_write_control_2c(filename,conf)
	elif conf['mode'] == 'SY':
		_write_control_mc(filename,conf,"syn")
	else:
		raise ValueError(f'[write_control] Unknown mode {mode["mode"]}')

def _write_control_1c(filename : str, conf : dict):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameters
	----------
	filename : str
		Save filename of the control file. Typically it is inv.trol
	config : dict
		Dictionary with the information from the config file
	
	

	"""
	import definitions as d
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
	
	if conf['psf'] == '':
		psf = ''
	else:
		psf = d.psf

	# Write lines
	with open(filename, 'w') as f:
		f.write(f'Number of cycles           (*):{cycles}                  ! (0=synthesis)\n')
		f.write('Observed profiles          (*):' + d.profile_obs + '      ! \n')
		f.write('Stray light file              :                   ! (none=no stray light contam)\n')
		f.write('PSF file                      :' + psf + '        ! (none=no convolution with PSF)\n')
		f.write('Wavelength grid file       (s):' + d.Grid + '! (none=automatic selection)\n')
		f.write('Atomic parameters file        :' + line + '    ! (none=DEFAULT LINES file)\n')
		f.write('Abundances file               :' + abundance + '         ! (none=DEFAULT ABUNDANCES file)\n')
		f.write('Initial guess model 1      (*):' + d.model_inv + '      !\n')
		f.write('Initial guess model 2         :\n')
		f.write('Weight for Stokes I           :' + weights[0] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes Q           :' + weights[1] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes U           :' + weights[2] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('Weight for Stokes V           :' + weights[3] + '                   ! (DEFAULT=1; 0=not inverted\n')
		f.write('AUTOMATIC SELECT. OF NODES?   :                   ! (DEFAULT=0=no; 1=yes)\n')
		f.write('Nodes for temperature 1       :' + nodes_temp + '\n')
		f.write('Nodes for electr. press. 1    :                         \n')
		f.write('Nodes for microturb. 1        :                         \n')
		f.write('Nodes for magnetic field 1    :' + nodes_magn + '\n')
		f.write('Nodes for LOS velocity 1      :' + nodes_vlos + '\n')
		f.write('Nodes for gamma 1             :' + nodes_gamma + '\n')
		f.write('Nodes for phi 1               :' + nodes_phi + '\n')
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
		f.write('mu=cos (theta)                :' + mu_cos +  '              ! (DEFAULT: mu=1)\n')
		f.write('Estimated S/N for I           :200                ! (DEFAULT: 1000) \n')
		f.write('Continuum contrast            :                    ! (DEFAULT: not used)\n')
		f.write('Tolerance for SVD             :' + d.SVD + '              ! (DEFAULT value: 1e-4)\n')
		f.write('Initial diagonal element      :                    ! (DEFAULT value: 1.e-3)\n')
		f.write('Splines/Linear Interpolation  :                    ! (0 or blank=splines, 1=linear)\n')
		f.write('Gas pressure at surface 1     :' + gas_pressure + '              ! (0 or blank=Pe boundary cond.)\n')
		f.write('Gas pressure at surface 2     :                    ! (0 or blank=Pe boundary cond.\n')
		f.write('Magnetic pressure term?       :                    ! (0 or blank=no, 1=yes\n')
		f.write("NLTE Departures filename      :                    ! blanck= LTE (Ej.) depart_6494.dat'\n")


	f.close()

def _write_control_2c(filename : str, conf : dict):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameters
	----------
	filename : str
		Save filename of the control file. Typically it is inv.trol
	conf : dict
		Dictionary with the information from the config file
	
	
	"""
	import definitions as d
	model1		= d.guess1			# Base Model 1
	model2		= d.guess2			# Base Model 2
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

	if " " in gas_pressure:
		gas_pressure1, gas_pressure2 = gas_pressure.split(" ")
	else:
		gas_pressure1  = gas_pressure2 = gas_pressure

	fill = conf["invert_fill"] # invert filling factor

	if conf['psf'] == '':
		psf = ''
	else:
		psf = d.psf

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
		f.write('Invert filling factor?        :' + fill  + '                    ! (0 or blank=no, 1=yes)\n')
		f.write('Invert stray light factor?    :0                   ! (0 or blank=no, 1=yes)\n')
		f.write('mu=cos (theta)                :'  + mu_cos +  '              ! (DEFAULT: mu=1)\n')
		f.write('Estimated S/N for I           :200                ! (DEFAULT: 1000) \n')
		f.write('Continuum contrast            :                    ! (DEFAULT: not used)\n')
		f.write('Tolerance for SVD             :' + d.SVD + '              ! (DEFAULT value: 1e-4)\n')
		f.write('Initial diagonal element      :                    ! (DEFAULT value: 1.e-3)\n')
		f.write('Splines/Linear Interpolation  :                    ! (0 or blank=splines, 1=linear)\n')
		f.write('Gas pressure at surface 1     :' + gas_pressure1 + '              ! (0 or blank=Pe boundary cond.)\n')
		f.write('Gas pressure at surface 2     :' + gas_pressure2 + '              ! (0 or blank=Pe boundary cond.\n')
		f.write('Magnetic pressure term?       :                    ! (0 or blank=no, 1=yes\n')
		f.write("NLTE Departures filename      :                    ! blanck= LTE (Ej.) depart_6494.dat'\n")


	f.close()


def _write_control_mc(filename : str, conf : dict, Type : str= 'inv'):
	"""
	Writes a control file in the scheme SIR expects it.
	
	Parameters
	----------
	filename : string
		Save filename of the control file. Typically it is inv.trol
	config : dict
		Dictionary with the information from the config file
	Type : string
		which type of control file is created ('syn' for synthesis, 'inv' for inversion)
	
	
	"""
	import definitions as d
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
		f.write('Tolerance for SVD             :' + d.SVD + '               ! (DEFAULT value: 1e-4)\n')
		f.write('Initial diagonal element      :                    ! (DEFAULT value: 1.e-3)\n')
		f.write('Splines/Linear Interpolation  :                    ! (0 or blank=splines, 1=linear)\n')
		f.write('Gas pressure at surface 1     :' + gas_pressure + '          ! (0 or blank=Pe boundary cond.)\n')
		f.write('Gas pressure at surface 2     :                    ! (0 or blank=Pe boundary cond.\n')
		f.write('Magnetic pressure term?       :                    ! (0 or blank=no, 1=yes\n')
		f.write("NLTE Departures filename      :                    ! blanck= LTE (Ej.) depart_6494.dat'\n")


	f.close()

def write_gauss_psf(sigma : float, filename : str) -> None:
	"""
	Writes the spectral point spread function with the given sigma. This function is used when the field psf in the config
	contains 'gauss 1.0' with sigma = 1.0 mA in this example.

	Parameters
	----------
	sigma : float
		Sigma of the Gauss function in mA
	filename : str
		Output file name

	

	"""
	Delta_ll = 20 # mA

	## Wavelength range to be printed in mA
	ll = np.arange(-(Delta_ll*31),(Delta_ll*32),Delta_ll)

	# Compute Gaussian
	g = np.exp(-(ll - 0)**2 / (2 * sigma**2))
	g = g / (sigma * np.sqrt(2 * np.pi)) # Normalised so that int g = 1

	f = open(filename, 'w')
	for num1,num2 in zip(ll,g):
		f.write(f" {num1:>9.4f}   {num2:>10.4E}\n")


def write_grid(conf : dict, filename : str = 'Grid.grid', waves : bool=None) -> None:
	"""
	Writes the Grid file with data from the config file

	Parameters
	----------
	config : dict
		Dictionary containing all the information from the config file
	filename : str,optional
		String containing the name of the Grid file, by default "Grid.grid"
	waves : numpy array,optional
		Array with the wavelength needed for mode '1C' and '2C', by default None

	Raises
	------
	ValueError
		if 'waves' not defined but needed for mode '1C' and '2C'
	ValueError
		if unknown mode

	
	"""
	if conf['mode'] == "MC" or conf['mode'] == "SY":
		_write_grid_mc(conf, filename)
	elif conf['mode'] == "1C":
		_write_grid(conf, filename, waves)
	elif conf['mode'] == "2C":
		_write_grid(conf, filename, waves)
	else:
		raise ValueError(f"[write_grid] Unknown Mode '{conf['mode']}'")

def _write_grid(conf : dict, filename : str, waves : np.array):
	"""
	Writes the Grid file with data from the config file

	Parameters
	----------
	config : dict
		Dictionary containing all the information from the config file
	filename : str
		String containing the name of the Grid file.
	waves : numpy.array
		Array with the wavelength

	Raises
	------
	ValueError
		if 'waves' not defined

	
	"""
	if waves is None:
		raise ValueError("[write_grid] 'waves' not defined")
	
	# Load data from config
	range_wave = conf['range_wave']
	line = read_line(os.path.join(conf['path'],conf['line']))
	atoms = conf["atoms"]

	# Define minimum, step and maximum
	Line_min = np.zeros(0, dtype=np.float64)
	Line_max = np.zeros(0, dtype=np.float64)
	Line_step = np.float64(conf["range_wave"][:,1]) # in mA
	
	# Correct if only one wavelength range is used to be an array
	if range_wave.shape[0] == 1:
		Line_step = np.array([Line_step])

	# Determine minimum and maximum to be used in sir
	for i in range(range_wave.shape[0]):
		Line_min  = np.append(Line_min,waves[np.argmin(np.abs(waves-range_wave[i,0]))])
		Line_max  = np.append(Line_max,Line_min[i] + Line_step[i]/1e3*(range_wave[i,2]-1))

	# Define wavelength grid to be saved
	with open(filename, 'w') as f:
		for i in range(len(atoms)):
			ind = np.where(line['Line'] == int(atoms[i].split(',')[0]))[0][0] # Which index in line file corresponds to the atom
			llambdas = (np.array([Line_min[i],Line_max[i]], dtype=np.float64) - line['wavelength'][ind])*1e3 # in mA ; Determine relative wavelengths
			f.write(f"{atoms[i]}: {'%6.4f' % llambdas[0]},     {'%2.6f' % Line_step[i]},     {'%6.4f' % llambdas[-1]}\n")
	

def _write_grid_mc(conf : dict, filename : str) -> None:
	"""
	Writes the Grid file with data from the config file for the MC or SY mode

	Parameters
	----------
	conf : dict
		Dictionary containing all the information from the config file
	filename : str, optional
		String containing the name of the Grid file, by default "Grid.grid"
	
	

	"""

	# Load data from config
	range_wave = conf['range_wave']
	atoms = conf["atoms"]


	# Define minimum, step and maximum
	Line_min  = np.array(range_wave[:,0]).astype(np.float64)
	Line_step = np.array(range_wave[:,1]).astype(np.float64)
	Line_max  = Line_min + Line_step*(np.array(range_wave[:,2]).astype(np.float64)-1)

	# Define wavelength grid to be saved
	with open(filename, 'w') as f:
		for i in range(len(atoms)):
			f.write(f"{atoms[i]}: {'%6.4f' % Line_min[i]},     {'%2.4f' % Line_step[i]},     {'%6.4f' % Line_max[i]}\n")

def write_model(filename, Header, log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z = None, Pg = None, rho = None):
	"""
	Write a model with the given data in a specific format. Note that negative values
	have one white space less

	Parameters
	----------
	filename : str
		Name of the saved file
	Header : str
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


def x_y_add_zeros(x : float, y : float) -> Tuple[str, str]:
	"""
	Adds zeros so that the returning strings have 4 letters

	Parameters
	----------
	x : float
		x position
	y : float
		y position

	Returns
	-------
	str
		x as a string of 4 letters
	str
		y as a string of 4 letters

	

	"""
	if x < 10:
		x_str = "000" + str(x)
	elif x < 100:
		x_str = "00" + str(x)
	elif x < 1000:
		x_str = "0" + str(x)
	else:
		x_str = str(x)
	if y < 10:
		y_str = "000" + str(y)
	elif y < 100:
		y_str = "00" + str(y)
	elif y < 1000:
		y_str = "0" + str(y)
	else:
		y_str = str(y)

	return x_str, y_str
