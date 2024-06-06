#!/usr/bin/env python3
"""
Function to create configuration files which can be executed directly from the terminal.
The function asks questions which needs to be answered and creates then a config file.
"""

import sys 
import sir
from os.path import exists

def help():
	print("create_config.py - Create a config file for the inversion")
	print("Usage: python create_config.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Config File")
	sir.option("-path","Set the path via terminal")
	print()
	sys.exit()

def _config_MC():
	"""
	Creates a config file for the Monte-Carlo simulation by asking questions.

	Parameters
	----------
	None

	Returns
	-------
	None
	"""
	File = input("File name: ")

	if exists(File):
		print('File exists already.')
		sys.exit()

	mode = "MC"
	if "-path" in sys.argv:
		path = sys.argv[sys.argv.index("-path")+1]
	else:
		path = input("Path: ")

	num	= input("Number of models: ")
	cycles = input("Cycles: ")
	model = input("Base model: ")
	model_nodes	= input("Create Models with 1, 2 or 3 nodes: ")
	
	model_pars = input("Randomize these parameters for creating models [B,T,vlos,gamma]: ")
	syn_out	= input("Synthesis prefix [syn]: ")
	noise_out = input("Noise Prefix [noise]: ")
	instrument = input("Instrument (GRIS, Hinode or blank) for predefined noise levels: ")
	if instrument == "GRIS":
		noise_I=0
		noise_Q=76.80e-5 # Computed by using many different quiet sun observations from GRIS
		noise_U=76.77e-5 # Computed by using many different quiet sun observations from GRIS
		noise_V=74.49e-5 # Computed by using many different quiet sun observations from GRIS
	elif instrument == "Hinode":
		noise_I=0
		noise_Q=1e-3
		noise_U=1e-3
		noise_V=1e-3
	else:
		noise_I = input("Noise in I: ")
		noise_Q = input("Noise in Q: ")
		noise_U = input("Noise in U: ")
		noise_V = input("Noise in V: ")
	inv_out		= input("Inversion output      [inv]: ")
	atoms		= input("Atoms (e.g. 8,9;3,4   ';' == newline): ")
	if instrument != 'Hinode':
		range_wave = input("Ranges in the wavelengths (in relative mA) to be considered (as 'min1,step1,max1;min2,step2,max2;...', ;=newline):" )
	else:
		range_wave = '-750.0,22.5,1750.0'
	random_guess = input("Number of random guess models (0 = use base model): ")
	if random_guess != '0':
		random_pars = input("Randomize these parameters [B,T,vlos,gamma]: ")
	else:
		random_pars = ''
	guess = input ("Take npy file as initial guess? Write name of the file, if used: ")

	weights = input ("Weights as a list (I,Q,U,V)   [1,1,1,1]: ")
	nodes_temp = input("Nodes in temperature  (as list)        : ")
	nodes_magn = input("Nodes in magn. field  (as list)        : ")
	nodes_vlos = input("Nodes in velocity los (as list)        : ")
	nodes_gamma = input("Nodes in inclination/gamma (as list)   : ")
	nodes_phi = input("Nodes in azimuth/phi (as list)         : ")
	vmacro		= input ("Value for the macroturbulence [0.1000]: ")
	abundance = input("Abundance file               [THEVENIN]: ")
	chi2 = input("Compute chi2 and save it under this file [chi2.bin]: ")
	line = input("Lines file                      [Lines]: ")
	gas_pressure = input("Gas Pressure Boundary condition  [-1 => 3.944e+3]: ")


	lim_B = input("Limits for randomising the magn. field in G             [0,4000]: ")
	lim_vlos = input("Limits for randomising the line-of-sight vel. in cm/s [-3e5,3e5]: ")
	lim_gamma = input("Limits for randomising the inclination in deg            [0,180]: ")
	lim_phi = input("Limits for randomising the azimuth in deg                [0,180]: ")

	if model_nodes == 1:
		print("With 1 node, the part after ';' can be ignored as it is not used in this setting.")
	# For 3 Nodes
	create_B = input("Limits for the magnetic field in G for creating models              [500,4000;0,1000]: ")
	create_vlos = input("Limits for the line-of-sight vel. in cm/s for creating models        [-3e5,3e5;-3e5,3e5]: ")
	create_gamma = input("Limits for gamma in deg for creating models                                [0,180;0,180]: ")
	create_phi = input("Limits for phi in deg for creating models                                  [0,180;0,180]: ")
	create_points = input("The nodes where the model values  are randomised and interpolated (increasing) [-4,-1,1]: ")

	if create_B == "":
		create_B = "500,4000;0,1000" # Note that there are add. conditions in the script handling strong fields
	if create_vlos == "":
		create_vlos = "-3e5,3e5;-3e5,3e5" # The 2nd part depends on the 1st and last point and is between those
	if create_gamma == "":
		create_gamma = "0,180;0,180" # The 2nd part depends on the 1st and last point and is between those
	if create_phi == "":
		create_phi = "0,180;0,180" # The 2nd part depends on the 1st and last point and is between those
	if create_points == "":
		create_points = "-4,-1,1" # At this log tau points the models are interpolated with splines (increasing)

	# Fill with standard values
	if lim_B == '':
		lim_B = '0,4000'
	if lim_vlos == '':
		lim_vlos = '-3e5,3e5'
	if lim_gamma == '':
		lim_gamma = '0,180'
	if lim_phi == '':
		lim_phi = '0,180'
	if syn_out == '':
		syn_out = 'syn'
	if noise_out == '':
		noise_out = 'noise'
	if inv_out == '':
		inv_out = 'inv'

	if abundance == '':
		abundance = 'THEVENIN'
	if vmacro == '':
		vmacro = '0.1000'
	if chi2 == '':
		chi2 = 'chi2.bin'
	if line == '':
		line = 'Lines'
	if weights == '':
		weights = '1,1,1,1'
	if random_pars == '':
		random_pars = "B,T,vlos,gamma"
	if model_pars == '':
		model_pars = "B,T,vlos,gamma"
	if gas_pressure == "-1":
		gas_pressure = "3.944e+3"



	conf = {
		"mode" : mode,
		"path" : path,
		"num" : num,
		"model" : model,
		"atoms" : atoms,
		"range_wave" : range_wave,
		"random_guess" : random_guess,
		"random_pars" : random_pars,
		"model_pars" : model_pars,
		"model_nodes" : model_nodes,
		"syn_out" : syn_out,
		"noise_out" : noise_out,
		"noise_I" : noise_I,
		"noise_Q" : noise_Q,
		"noise_U" : noise_U,
		"noise_V" : noise_V,
		"create_B" : create_B,
		"create_vlos" : create_vlos,
		"create_gamma" : create_gamma,
		"create_phi" : create_phi,
		"create_points" : create_points,
		"inv_out" : inv_out,
		"chi2" : chi2,
		"line" : line,
		"guess" : guess,
		"cycles" : cycles,
		"weights" : weights,
		"nodes_temp" : nodes_temp,
		"nodes_magn" : nodes_magn,
		"nodes_vlos" : nodes_vlos,
		"nodes_gamma" : nodes_gamma,
		"nodes_phi" : nodes_phi,
		"vmacro" : vmacro,
		"abundance" : abundance,
		"gas_pressure" : gas_pressure,
		"lim_B" : lim_B,
		"lim_vlos" : lim_vlos,
		"lim_gamma" : lim_gamma,
		"lim_phi" : lim_phi,
	}

	sir.write_config(File,conf)

	return

def _config_1C():
	"""
	Creates a config file for the 1 Component Inversion by asking questions.

	Parameters
	----------
	None

	Returns
	-------
	None
	"""
	File = input("File name: ")

	if exists(File):
		print('File exists already.')
		sys.exit()

	mode = "1C"
	if "-path" in sys.argv:
		path = sys.argv[sys.argv.index("-path")+1]
	else:
		path		= input ("Path, where the files are: ")
	#cube			= input ("Location of the Data cube for preprocessing (format is nx,ny,ns,nwave) as a .bin file (can be left empty if no preprocessing): ")
	cube_inv		= input ("Location of the Data cube used for the inversion (format is nx,ny,ns,nwave) in the path as a bin file: ")
	preprocess		= input ("Preprocess data? (Normalisation and/or spectral veil correction? (yes -> 1, no -> 0, perform directly inversion): ")
	Map			= input ("Map in pixels to be used for the inversion (format xmin,xmax,ymin,ymax): ")

	instrument	= input ("Instrument          (GRIS, Hinode or blank): ")
	# Ask for spectral veil for Gris:
	if instrument == 'Hinode':
		fts_file = ''
	else:
		fts_file	= input ("Absolute path to fts file (if not blank => correct spectral veil): ")
	shift_wave	= input ("Shift the wavelength grid in mA: [0]: ")
	save_cube	= input ("Save preprocessed data? (1 = True, 0 = False): ")
	
	quiet_sun		= input ("Quiet sun region as a list (format x1,x2,y1,y2; 0 = already normalised): ")
	cycles		= input ("Cycles: ")
	model		= input ("Base model: ")
	inv_out		= input ("Inversion Prefix [out]: ")
	line			= input ("Line file                       [Lines]: ")
	atoms			= input ("Atoms (e.g. 8,9;3,4   ';' == newline): ")
	range_wave		= input ("Range for the grid file as (Start wavelength in abs. wavelength, Step in mA, Number of wavelenghts) for each line in the grid file:" )
	random_guess	= input ("Number of random guess models (0 = use base model): ")
	if random_guess != '0':
		random_pars	= input ("Randomize these parameters [B,T,vlos,gamma]: ")
	else:
		random_pars = ''
	guess		= input ("Take bin file as initial guess? Write name of the file, if used: ")
	psf			 = input ("Filename of psf (.dat file) or 'gauss xx.xx' with sigma=xx.xx in mA, blank = not used): ")

	weights		= input ("Weights as a list (I,Q,U,V)   [1,1,1,1]: ")
	nodes_temp	= input ("Nodes in temperature  (as list)        : ")
	nodes_magn	= input ("Nodes in magn. field  (as list)        : ")
	nodes_vlos	= input ("Nodes in velocity los (as list)        : ")
	nodes_gamma	= input ("Nodes in inclination/gamma (as list)   : ")
	nodes_phi		= input ("Nodes in azimuth/phi (as list)    	 : ")
	vmacro		= input ("Value for the macroturbulence [0.1000]: ")
	mu_cos		= input ("mu = cos theta                         : ")
	abundance		= input ("Abundance file               [THEVENIN]: ")
	chi2 = input("Compute chi2 and save it under this file [out_chi2.bin]: ")
	gas_pressure   = input ("Gas Pressure Boundary condition  [-1 => 3944]: ")

	lim_B		= input ("Limits for randomising the magn. field in G             [0,5000]: ")
	lim_vlos		= input ("Limits for randomising the line-of-sight vel. in cm/s [-1e5,1e5]: ")
	lim_gamma		= input ("Limits for randomising the inclination in deg            [0,180]: ")
	lim_azimuth	= input ("Limits for randomising the azimuth in deg                [0,180]: ")



	# Fill with standard values
	if shift_wave == '':
		shift_wave = '0'
	if lim_B == '':
		lim_B = '0,5000'
	if lim_vlos == '':
		lim_vlos = '-1e5,1e5'
	if lim_gamma == '':
		lim_gamma = '0,180'
	if lim_azimuth == '':
		lim_azimuth = '0,180'
	if inv_out == '':
		inv_out = 'out'
	if vmacro == '':
		vmacro = "0.1000"
	if line == '':
		line = 'Lines'
	if abundance == '':
		abundance = 'THEVENIN'
	if chi2 == '':
		chi2 = 'out_chi2.bin'
	if weights == '':
		weights = '1,1,1,1'
	if random_pars == '':
		random_pars = "B,T,vlos,gamma"
	if gas_pressure == "-1":
		gas_pressure = "3.944"

	conf = {
		"mode" : mode,
		"path" : path,
		#"cube" : cube,
		"cube_inv" : cube_inv,
		"map" : Map,
		"instrument" : instrument,
		"shift_wave" : shift_wave,
		"save_cube" : save_cube,
		"preprocess" : preprocess,
		"quiet_sun" : quiet_sun,
		"fts_file" : fts_file,
		"model" : model,
		"range_wave" : range_wave,
		"inv_out" : inv_out,
		"chi2" : chi2,
		"line" : line,
		"atoms" : atoms,
		"guess" : guess,
		"psf" : psf,
		"cycles" : cycles,
		"weights" : weights,
		"nodes_temp" : nodes_temp,
		"nodes_magn" : nodes_magn,
		"nodes_vlos" : nodes_vlos,
		"nodes_gamma" : nodes_gamma,
		"nodes_phi" : nodes_phi,
		"vmacro" : vmacro,
		"mu_cos" : mu_cos,
		"abundance" : abundance,
		"gas_pressure" : gas_pressure,
		"random_guess" : random_guess,
		"random_pars" : random_pars,
		"lim_B" : lim_B,
		"lim_vlos" : lim_vlos,
		"lim_gamma" : lim_gamma,
		"lim_phi" : lim_azimuth,
	}
	sir.write_config(File,conf)

def _config_2C():
	"""
	Creates a config file for the 2 Components Inversion by asking questions.

	Parameters
	----------
	None

	Returns
	-------
	None
	"""
	File = input("File name: ")

	if exists(File):
		print('File exists already.')
		sys.exit()

	mode = "2C"
	if "-path" in sys.argv:
		path = sys.argv[sys.argv.index("-path")+1]
	else:
		path		= input ("Path: ")
	#cube			= input ("Location of the Data cube for preprocessing (format is nx,ny,ns,nwave) as a fits or npy file (can be left empty if no preprocessing): ")
	cube_inv		= input ("Location of the Data cube used for the inversion (format is nx,ny,ns,nwave) in the path as a fits or npy file: ")
	preprocess		= input ("Preprocess data? (Normalisation and/or spectral veil correction? (yes -> 1, no -> 0, perform directly inversion): ")
	Map			= input ("Map in pixels (format xmin,xmax,ymin,ymax, 0 => all pixels): ")

	instrument	= input ("Instrument          (GRIS, Hinode or blank): ")
	# Ask for spectral veil for Gris:
	if instrument == 'Hinode':
		fts_file = ''
	else:
		fts_file	= input ("Absolute fts file (if not blank => correct spectral veil): ")
	shift_wave	= input ("Shift the wavelength grid in mA: [0]: ")
	save_cube	= input ("Save preprocessed data? (1 = True, 0 = False): ")

	quiet_sun		= input ("Quiet sun region as a list (format x1,x2,y1,y2; 0 = already normalised): ")
	cycles		= input ("Cycles: ")
	model1		= input ("Base model 1: ")
	model2		= input ("Base model 2: ")
	inv_out		= input ("Inversion output prefix [out]: ")
	line			= input ("Line file                       [Lines]: ")
	atoms			= input ("Atoms (e.g. 8,9;3,4   ';' == newline): ")
	range_wave  = input("Range for the grid file as (Start wavelength in abs. wavelength, Step in mA, Number of wavelenghts) for each line in the grid file:")
	random_guess	= input ("Number of random guess models (0 = use base model): ")
	if random_guess != '0':
		random_pars	= input ("Randomize these parameters [B,T,vlos,gamma]: ")
		fill		= input ("Filling factor for the two models? Seperated by a ',': ")
	else:
		random_pars = ''
		fill = ''
	guess1		= input ("Take bin file as initial guess for model 1? Write name of the file, if used: ")
	guess2		= input ("Take bin file as initial guess for model 2? Write name of the file, if used: ")
	psf			 = input ("Filename of psf (.dat file) or 'gauss xx.xx' with sigma=xx.xx in mA, blank = not used: ")
	fill		= input ("Filling factor for the two models? Seperated by a ',': ")

	weights		= input ("Weights as a list (I,Q,U,V)   [1,1,1,1]: ")
	nodes_temp1	= input ("Nodes 1 in temperature  (as list)		: ")
	nodes_magn1	= input ("Nodes 1 in magn. field  (as list)        : ")
	nodes_vlos1	= input ("Nodes 1 in velocity los (as list)        : ")
	nodes_gamma1	= input ("Nodes 1 in inclination/gamma (as list)   : ")
	nodes_phi1	= input ("Nodes 1 in azimuth/phi (as list)         : ")
	nodes_temp2	= input ("Nodes 2 in temperature  (as list)        : ")
	nodes_magn2	= input ("Nodes 2 in magn. field  (as list)        : ")
	nodes_vlos2	= input ("Nodes 2 in velocity los (as list)        : ")
	nodes_gamma2	= input ("Nodes 2 in inclination/gamma (as list)   : ")
	nodes_phi2	= input ("Nodes 2 in azimuth/phi (as list)         : ")
	vmacro		= input ("Value for the macroturbulence [0.1000]: ")
	mu_cos		= input ("mu = cos theta                         : ")
	abundance		= input ("Abundance file               [THEVENIN]: ")
	chi2 = input("Compute chi2 and save it under this file [out_chi2.bin]: ")
	gas_pressure   = input ("Gas Pressure Boundary condition  [-1 => 3.944e+3]: ")

	lim_B1		= input ("Limits 1 for randomising the magn. field in G             [0,5000]: ")
	lim_vlos1		= input ("Limits 1 for randomising the line-of-sight vel. in cm/s [-1e5,1e5]: ")
	lim_gamma1	= input ("Limits 1 for randomising the inclination in deg            [0,180]: ")
	lim_azimuth1	= input ("Limits 1 for randomising the azimuth in deg                [0,180]: ")
	lim_B2		= input ("Limits 2 for randomising the magn. field in G             [0,5000]: ")
	lim_vlos2		= input ("Limits 2 for randomising the line-of-sight vel. in cm/s [-1e5,1e5]: ")
	lim_gamma2	= input ("Limits 2 for randomising the inclination in deg            [0,180]: ")
	lim_azimuth2	= input ("Limits 2 for randomising the azimuth in deg                [0,180]: ")




	# Fill with standard values
	if shift_wave == '':
		shift_wave = '0'
	if lim_B1 == '':
		lim_B1 = '0,5000'
	if lim_vlos1 == '':
		lim_vlos1 = '-1e5,1e5'
	if lim_gamma1 == '':
		lim_gamma1 = '0,180'
	if lim_azimuth1 == '':
		lim_azimuth1 = '0,180'
	if lim_B2 == '':
		lim_B2 = '0,5000'
	if lim_vlos2 == '':
		lim_vlos2 = '-1e5,1e5'
	if lim_gamma2 == '':
		lim_gamma2 = '0,180'
	if lim_azimuth2 == '':
		lim_azimuth2 = '0,180'

	if inv_out == '':
		inv_out = 'out'
	if line == '':
		line = 'Lines'
	if abundance == '':
		abundance = 'THEVENIN'
	if chi2 == '':
		chi2 = 'out_chi2.bin'
	if weights == '':
		weights = '1,1,1,1'
	if vmacro == '':
		vmacro = '0.1000'
	if random_pars == '':
		random_pars = "B,T,vlos,gamma"
	if gas_pressure == "-1":
		gas_pressure = "3.944e+3"

	conf = {
		"mode" : mode,
		"path" : path,
		#"cube" : cube,
		"cube_inv" : cube_inv,
		"map" : Map,
		"instrument" : instrument,
		"shift_wave" : shift_wave,
		"preprocess" : preprocess,
		"quiet_sun" : quiet_sun,
		"fts_file" : fts_file,
		"save_cube" : save_cube,
		"model1" : model1,
		"model2" : model2,
		"fill" : fill,
		"range_wave" : range_wave,
		"inv_out" : inv_out,
		"chi2" : chi2,
		"line" : line,
		"atoms" : atoms,
		"guess1" : guess1,
		"guess2" : guess2,
		"psf" : psf,
		"cycles" : cycles,
		"weights" : weights,
		"nodes_temp1" : nodes_temp1,
		"nodes_magn1" : nodes_magn1,
		"nodes_vlos1" : nodes_vlos1,
		"nodes_gamma1" : nodes_gamma1,
		"nodes_phi1" : nodes_phi1,
		"nodes_temp2" : nodes_temp2,
		"nodes_magn2" : nodes_magn2,
		"nodes_vlos2" : nodes_vlos2,
		"nodes_gamma2" : nodes_gamma2,
		"nodes_phi2" : nodes_phi2,
		"vmacro" : vmacro,
		"mu_cos" : mu_cos,
		"abundance" : abundance,
		"gas_pressure" : gas_pressure,
		"random_guess" : random_guess,
		"random_pars" : random_pars,
		"lim_B1" : lim_B1,
		"lim_vlos1" : lim_vlos1,
		"lim_gamma1" : lim_gamma1,
		"lim_azimuth1" : lim_azimuth1,
		"lim_B2" : lim_B2,
		"lim_vlos2" : lim_vlos2,
		"lim_gamma2" : lim_gamma2,
		"lim_azimuth2" : lim_azimuth2,
	}

	sir.write_config(File, conf)

def create_config():
	"""
	Creates a configuration file for the different modes
	"""
	mode = input ("Which mode do you want to use? [1C/2C/MC]: ")
	if mode == "1C":
		_config_1C()
	elif mode == "2C":
		_config_2C()
	elif mode == "MC":
		_config_MC()
	else:
		print("Mode unknown")

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()

	create_config()

