"""
Create a config file as expected for the inversion
"""

import numpy as np 
import sys 
import sir
from os.path import exists
import os

def help():
	print("create_config.py - Create a config file for the inversion")
	print("Usage: python create_config.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Config File")
	sir.option("-path","Set the path via terminal")
	print()
	sys.exit()

def config_MC():
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
	model_nodes	= input("Create Models with 2 or 3 nodes: ")
	model_out = input("Name of the npy file with the models [models_real.npy]: ")
	model_pars = input("Randomize these parameters for creating models [B,T,vlos,gamma]: ")
	syn_out	= input("Synthesis output [synthesis_profiles.npy]: ")
	noise_out = input("Noise output         [noise_profiles.npy]: ")
	instrument = input("Instrument (GRIS, Hinode or blank): ")
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
		noise_B = input("Noise in V: ")
	inv_out		= input("Inversion output      [inversion]: ")
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
	chi2 = input("File for chi2 map            [chi2.npy]: ")
	line = input("Lines file                      [Lines]: ")
	gas_pressure = input("Gas Pressure Boundary condition  [-1 => 3.944e+3]: ")


	lim_B = input("Limits for randomising the magn. field in G             [0,5000]: ")
	lim_vlos = input("Limits for randomising the line-of-sight vel. in cm/s [-3e5,3e5]: ")
	lim_gamma = input("Limits for randomising the inclination in deg            [0,180]: ")
	lim_phi = input("Limits for randomising the azimuth in deg                [0,180]: ")

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
		lim_B = '0,5000'
	if lim_vlos == '':
		lim_vlos = '-3e5,3e5'
	if lim_gamma == '':
		lim_gamma = '0,180'
	if lim_phi == '':
		lim_phi = '0,180'
	if model_out == '':
		model_out = 'models_real.npy'
	if syn_out == '':
		syn_out = 'synthesis_profiles.npy'
	if noise_out == '':
		noise_out = 'noise_profiles.npy'
	if inv_out == '':
		inv_out = 'inversion'

	if abundance == '':
		abundance = 'THEVENIN'
	if vmacro == '':
		abundance = '0.1000'
	if chi2 == '':
		chi2 = 'chi2.npy'
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



	with open(File, 'w') as f:
		f.write("# This is the config file, generated with create_config.py\n")
		f.write(f"mode : {mode} # Mode which code is executed\n")
		f.write(f"path : {path} # Path location where all the data is stored and will be saved\n")
		f.write(f"# \n")
		f.write(f"# General Stuff\n")
		f.write(f"# \n")
		f.write(f"num : {num} # Number of Models\n")
		f.write(f"instrument : {instrument} # Instrument used (GRIS, Hinode or empty)\n")
		f.write(f"model : {model} # Base Model for guess\n")
		f.write(f"atoms : {atoms} # Atoms to be used in Grid file\n")	
		f.write(f"range_wave : {range_wave} # Ranges of wavelengths in mA to be considered min1,step,max1;min2,step,max2;... First pair belongs to first line in Grid file, etc.\n")
		f.write(f"random_guess : {random_guess} # Create random guesses, 0 = use model as guess\n")
		f.write(f"random_pars : {random_pars} # Randomise these parameters for the guess\n")
		f.write(f"model_pars : {model_pars} # Randomize these parameters as a list\n")

		f.write(f"#\n")
		f.write(f"# Creating Models and Synthesis\n")
		f.write(f"#\n")
		f.write(f"model_nodes : {model_nodes} # Create models with 2 or 3 nodes\n")
		f.write(f"model_out : {model_out} # Output file of the created models as npy\n")
		f.write(f"syn_out : {syn_out} # Output of the synthesis profiles and models\n")
		f.write(f"noise_out : {noise_out} # Output of the noise profiles as npy\n")
		f.write(f"noise_I : {noise_I} # Noise in Q\n")	
		f.write(f"noise_Q : {noise_Q} # Noise in Q\n")
		f.write(f"noise_U : {noise_U} # Noise in U\n")
		f.write(f"noise_V : {noise_V} # Noise in V\n")
		f.write(f"create_B : {create_B} # The limits for the first and last node in B\n")
		f.write(f"create_vlos : {create_vlos} # The limits for the first and last node in vlos\n")
		f.write(f"create_gamma : {create_gamma} # The limits for the first and last node in gamma\n")
		f.write(f"create_phi : {create_phi} # The limits for the first and last node in phi\n")
		f.write(f"create_points : {create_points} # At this log tau points the models are interpolated with splines (increasing), 2 or 3 values for 2 or 3 nodes\n")

		f.write(f"# \n")
		f.write(f"# Inversion configuration\n")
		f.write(f"# \n")
		f.write(f"inv_out : {inv_out} # Prefix of the output of the inversion files\n")
		f.write(f"chi2 : {chi2} # Output of the chi2 values (npy)\n")
		f.write(f"line : {line} # Line file\n")
		f.write(f"guess : {guess} # Use a npy file as initial guesses, blank use base model\n")

		f.write(f"# \n")
		f.write(f"# Control file\n")
		f.write(f"# \n")
		f.write(f"cycles : {cycles} # Number of cycles\n")
		f.write(f"weights : {weights} # Weights in the control file\n")
		f.write(f"nodes_temp : {nodes_temp} # Nodes in T\n")
		f.write(f"nodes_magn : {nodes_magn} # Nodes in B\n")
		f.write(f"nodes_vlos : {nodes_vlos} # Nodes in vlos\n")
		f.write(f"nodes_gamma : {nodes_gamma} # Nodes in gamma\n")
		f.write(f"nodes_phi : {nodes_phi} # Nodes in phi\n")
		f.write(f"vmacro : {vmacro} # Macroturbulence velocity\n")
		f.write(f"abundance : {abundance} # Abundance file\n")
		f.write(f"gas_pressure : {gas_pressure} # Gas Pressure Boundary condition\n")

		f.write(f"# \n")
		f.write(f"# Radomisation Settings\n")
		f.write(f"# \n")
		f.write(f"lim_B : {lim_B} # Limits for the randomisation in B in G\n")
		f.write(f"lim_vlos : {lim_vlos} # Limits for the randomisation in vlos in cm/s\n")
		f.write(f"lim_gamma : {lim_gamma} # Limits for the randomisation in the inclination in deg\n")
		f.write(f"lim_phi : {lim_phi} # Limits for the randomisation in the azimuth in deg")
	

def config_1C():
    File = input("File name: ")

    if exists(File):
        print('File exists already.')
        sys.exit()

    mode = "1C"
    if "-path" in sys.argv:
        path = sys.argv[sys.argv.index("-path")+1]
    else:
        path		= input ("Path, where the files are: ")
    cube			= input ("Location of the Data cube for preprocessing (format is nx,ny,ns,nwave) as a npy file (can be left empty if no preprocessing): ")
    cube_inv		= input ("Location of the Data cube used for the inversion (format is nx,ny,ns,nwave) in the path as a bin file: ")
    preprocess		= input ("Preprocess data? (Normalisation and/or spectral veil correction? (yes -> 1, no -> 0, perform directly inversion): ")
    Map			= input ("Map in pixels (format xmin,xmax,ymin,ymax): ")

    instrument	= input ("Instrument          (GRIS, Hinode or blank): ")
    # Ask for spectral veil for Gris:
    if instrument == 'Hinode':
        fts_file = ''
    else:
        fts_file	= input ("Absolute fts file (if not blank => correct spectral veil): ")
    shift_wave	= input ("Shift the wavelength grid in mA: [0]: ")

    quiet_sun		= input ("Quiet sun region as a list (format x1,x2,y1,y2; 0 = already normalised): ")
    cycles		= input ("Cycles: ")
    model		= input ("Base model: ")
    inv_out		= input ("Inversion output as npy [inversion.npy]: ")
    line			= input ("Line file                       [Lines]: ")
    atoms			= input ("Atoms (e.g. 8,9;3,4   ';' == newline): ")
    Range		= input ("Ranges in the wavelengths (pixel or angstrom) to be considered (as 'min1,max1;min2,max2;...') defining the indices linewise as in grid file):" )
    Step		= input ("Wavelength steps in mA (as 'step1,step2,...') defining the indices linewise as in grid file):" )
    random_guess	= input ("Number of random guess models (0 = use base model): ")
    if random_guess != '0':
        random_pars    = input ("Randomize these parameters [B,T,vlos,gamma]: ")
    else:
        random_pars = ''
    guess		= input ("Take npy file as initial guess? Write name of the file, if used: ")
    psf		     = input ("Filename of psf (.dat file, if it does not exist => Compute from spectral veil corr., blank = not used): ")

    weights		= input ("Weights as a list (I,Q,U,V)   [1,1,1,1]: ")
    nodes_temp	= input ("Nodes in temperature  (as list)        : ")
    nodes_magn	= input ("Nodes in magn. field  (as list)        : ")
    nodes_vlos	= input ("Nodes in velocity los (as list)        : ")
    nodes_gamma	= input ("Nodes in inclination/gamma (as list)   : ")
    nodes_phi		= input ("Nodes in azimuth/phi (as list)         : ")
    vmacro		= input ("Value for the macroturbulence [0.1000]: ")
    mu_cos		= input ("mu = cos theta                         : ")
    abundance		= input ("Abundance file               [THEVENIN]: ")
    chi2			= input ("File for chi2 map            [chi2.npy]: ")
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
        inv_out = 'inversion.npy'
    if vmacro == '':
        vmacro = "0.1000"
    if line == '':
        line = 'Lines'
    if abundance == '':
        abundance = 'THEVENIN'
    if chi2 == '':
        chi2 = 'chi2.npy'
    if weights == '':
        weights = '1,1,1,1'
    if random_pars == '':
        random_pars = "B,T,vlos,gamma"
    if gas_pressure == "-1":
        gas_pressure = "3.944"

    with open(File, 'w') as f:
        f.write("# This is the config file, generated with create_config.py\n")
        f.write(f"mode : {mode} # Mode which code is executed\n")
        f.write(f"path : {path} # Path location where all the data is stored and will be saved\n")
        f.write(f"# \n")
        f.write(f"# Stuff from the data\n")
        f.write(f"# \n")
        f.write(f"cube : {cube} # Data cube name (npy or fits) used for preprocessing data if 'preprocess' is 1\n")
        f.write(f"cube_inv : {cube_inv} # Data cube name for the inversion (npy or fits)\n")
        f.write(f"map : {Map} # Pixels to be considered as a list\n")
        f.write(f"instrument : {instrument} # Instrument used (GRIS, Hinode or empty)\n")
        f.write(f"shift_wave : {shift_wave} # Shift the wavelength grid when waves file is created in mA\n")

        f.write("#\n")
        f.write("# Data Preprocessing when main.py is executed\n")
        f.write("# If the script is executed directly, it is not affected. Not that 'cube_inv' will be overwritten!\n")
        f.write("#\n")
        f.write(f"preprocess : {preprocess} # Preprocess data (1 = True, 0 = False)\n")
        f.write(f"quiet_sun : {quiet_sun} # Quiet sun region for normalization as a list (0 => already normalised)\n")
        f.write(f"fts_file : {fts_file} # FTS file, blank = do not correct spectral veil\n")

        f.write(f"# \n")
        f.write(f"# Inversion configuration\n")
        f.write(f"# \n")
        f.write(f"model : {model} # Base Model for guess\n")
        f.write(f"range_wave : {Range} # Ranges of wavelengths (pixel or Angstrom) to be considered min1,max1;min2,max2;... First pair belongs to first line in Grid file, etc.\n")
        f.write(f"step_wave : {Step} # Step between wavelength points in mA to be considered as Step1,Step2,... First value belongs to first line in Grid file, etc.\n")
        f.write(f"inv_out : {inv_out} # Prefix of the output of the inversion files\n")
        f.write(f"chi2 : {chi2} # Output of the chi2 values (npy)\n")
        f.write(f"line : {line} # Line file\n")
        f.write(f"atoms : {atoms} # Atoms to be used in Grid file\n")
        f.write(f"guess : {guess} # Use a npy file as initial guesses, blank use base model\n")
        f.write(f"psf : {psf} # .dat file (if it does not exist, computed from spectral veil parameter), blank=not used\n")
        f.write(f"# \n")
        f.write(f"# Control file\n")
        f.write(f"# \n")
        f.write(f"cycles : {cycles} # Number of cycles\n")
        f.write(f"weights : {weights} # Weights in the control file\n")
        f.write(f"nodes_temp : {nodes_temp} # Nodes in T\n")
        f.write(f"nodes_magn : {nodes_magn} # Nodes in B\n")
        f.write(f"nodes_vlos : {nodes_vlos} # Nodes in vlos\n")
        f.write(f"nodes_gamma : {nodes_gamma} # Nodes in gamma\n")
        f.write(f"nodes_phi : {nodes_phi} # Nodes in phi\n")
        f.write(f"vmacro : {vmacro} # Macroturbulence velocity\n")
        f.write(f"mu_cos : {mu_cos} # mu = cos theta\n")
        f.write(f"abundance : {abundance} # Abundance file\n")
        f.write(f"gas_pressure : {gas_pressure} # Gas Pressure Boundary condition\n")
        f.write(f"# \n")
        f.write(f"# Radomisation Settings\n")
        f.write(f"# \n")
        f.write(f"random_guess : {random_guess} # Create random guesses, 0 = use model as guess\n")
        f.write(f"random_pars : {random_pars} # Randomise these parameters in the file(s) below\n")
        f.write(f"lim_B : {lim_B} # Limits for the randomisation in B in G\n")
        f.write(f"lim_vlos : {lim_vlos} # Limits for the randomisation in vlos in cm/s\n")
        f.write(f"lim_gamma : {lim_gamma} # Limits for the randomisation in the inclination in deg\n")
        f.write(f"lim_phi : {lim_azimuth} # Limits for the randomisation in the azimuth in deg")
        

def config_2C():
    File = input("File name: ")

    if exists(File):
        print('File exists already.')
        sys.exit()

    mode = "2C"
    if "-path" in sys.argv:
        path = sys.argv[sys.argv.index("-path")+1]
    else:
        path		= input ("Path: ")
    cube			= input ("Location of the Data cube for preprocessing (format is nx,ny,ns,nwave) as a fits or npy file (can be left empty if no preprocessing): ")
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

    quiet_sun		= input ("Quiet sun region as a list (format x1,x2,y1,y2; 0 = already normalised): ")
    cycles		= input ("Cycles: ")
    model1		= input ("Base model 1: ")
    model2		= input ("Base model 2: ")
    inv_out		= input ("Inversion output prefix [inversion]: ")
    line			= input ("Line file                       [Lines]: ")
    atoms			= input ("Atoms (e.g. 8,9;3,4   ';' == newline): ")
    Range		= input ("Ranges in the wavelengths (pixel or angstrom) to be considered (as 'min1,max1;min2,max2;...') defining the indices linewise as in grid file):" )
    Step		= input ("Wavelength steps in mA (as 'step1,step2,...') defining the indices linewise as in grid file):" )
    random_guess	= input ("Number of random guess models (0 = use base model): ")
    if random_guess != '0':
        random_pars    = input ("Randomize these parameters [B,T,vlos,gamma]: ")
    else:
        random_pars = ''
    guess1		= input ("Take bin file as initial guess for model 1? Write name of the file, if used: ")
    guess2		= input ("Take bin file as initial guess for model 2? Write name of the file, if used: ")
    psf		     = input ("Filename of psf (.dat file, if it does not exist => Compute from spectral veil corr., blank = not used): ")

    weights		= input ("Weights as a list (I,Q,U,V)   [1,1,1,1]: ")
    nodes_temp1	= input ("Nodes 1 in temperature  (as list)        : ")
    nodes_magn1	= input ("Nodes 1 in magn. field  (as list)        : ")
    nodes_vlos1	= input ("Nodes 1 in velocity los (as list)        : ")
    nodes_gamma1	= input ("Nodes 1 in inclination/gamma (as list)   : ")
    nodes_phi1	= input ("Nodes 1 in azimuth/phi (as list)         : ")
    nodes_temp2	= input ("Nodes 2 in temperature  (as list)        : ")
    nodes_magn2	= input ("Nodes 2 in magn. field  (as list)        : ")
    nodes_vlos2	= input ("Nodes 2 in velocity los (as list)        : ")
    nodes_gamma2	= input ("Nodes 2 in inclination/gamma (as list)   : ")
    nodes_phi2	= input ("Nodes 2 in azimuth/phi (as list)         : ")

    mu_cos		= input ("mu = cos theta                         : ")
    abundance		= input ("Abundance file               [THEVENIN]: ")
    chi2			= input ("File for chi2 map            [chi2.npy]: ")
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
        inv_out = 'inversion'
    if line == '':
        line = 'Lines'
    if abundance == '':
        abundance = 'THEVENIN'
    if chi2 == '':
        chi2 = 'chi2.npy'
    if weights == '':
        weights = '1,1,1,1'
    if random_pars == '':
        random_pars = "B,T,vlos,gamma"
    if gas_pressure == "-1":
        gas_pressure = "3.944e+3"

    with open(File, 'w') as f:
        f.write("# This is the config file, generated with create_config.py\n")
        f.write(f"mode : {mode} # Mode which code is executed\n")
        f.write(f"path : {path} # Path location where all the data is stored and will be saved\n")
        f.write(f"# \n")
        f.write(f"# Stuff from the data\n")
        f.write(f"# \n")
        f.write(f"cube : {cube} # Data cube name (npy or fits) used for preprocessing data if 'preprocess' is 1\n")
        f.write(f"cube_inv : {cube_inv} # Data cube name for the inversion (npy or fits)\n")
        f.write(f"map : {Map} # Pixels to be considered as a list (0 means all pixels)\n")
        f.write(f"instrument : {instrument} # Instrument used (GRIS, Hinode or empty)\n")
        f.write(f"shift_wave : {shift_wave} # Shift the wavelength grid when waves file is created in mA\n")

        f.write("#\n")
        f.write("# Data Preprocessing when main.py is executed\n")
        f.write("# If the script is executed directly, it is not affected. Not that 'cube_inv' will be overwritten!\n")
        f.write("#\n")
        f.write(f"preprocess : {preprocess} # Preprocess data (1 = True, 0 = False)\n")
        f.write(f"quiet_sun : {quiet_sun} # Quiet sun region for normalization as a list (0 => already normalised)\n")
        f.write(f"fts_file : {fts_file} # FTS file, blank = do not correct spectral veil\n")

        f.write(f"# \n")
        f.write(f"# Inversion configuration\n")
        f.write(f"# \n")
        f.write(f"model1 : {model1} # Base Model 1 for guess\n")
        f.write(f"model2 : {model2} # Base Model 2 for guess\n")
        f.write(f"range_wave : {Range} # Ranges of wavelengths (pixel or Angstrom) to be considered min1,max1;min2,max2;... First pair belongs to first line in Grid file, etc.\n")
        f.write(f"step_wave : {Step} # Step between wavelength points in mA to be considered as Step1,Step2,... First value belongs to first line in Grid file, etc.\n")
        f.write(f"inv_out : {inv_out} # Prefix of the output of the inversion files\n")
        f.write(f"chi2 : {chi2} # Output of the chi2 values (npy)\n")
        f.write(f"line : {line} # Line file\n")
        f.write(f"atoms : {atoms} # Atoms to be used in Grid file\n")
        f.write(f"guess1 : {guess1} # Use a npy file as initial guesses, blank use base model 1\n")
        f.write(f"guess2 : {guess2} # Use a npy file as initial guesses, blank use base model 2\n")
        f.write(f"psf : {psf} # .dat file (if it does not exist, computed from spectral veil parameter), blank=not used\n")
        f.write(f"# \n")
        f.write(f"# Control file\n")
        f.write(f"# \n")
        f.write(f"cycles : {cycles} # Number of cycles\n")
        f.write(f"weights : {weights} # Weights in the control file\n")
        f.write(f"nodes_temp1 : {nodes_temp1} # Nodes in T\n")
        f.write(f"nodes_magn1 : {nodes_magn1} # Nodes in B\n")
        f.write(f"nodes_vlos1 : {nodes_vlos1} # Nodes in vlos\n")
        f.write(f"nodes_gamma1 : {nodes_gamma1} # Nodes in gamma\n")
        f.write(f"nodes_phi1 : {nodes_phi1} # Nodes in phi\n")
        f.write(f"nodes_temp2 : {nodes_temp2} # Nodes in T\n")
        f.write(f"nodes_magn2 : {nodes_magn2} # Nodes in B\n")
        f.write(f"nodes_vlos2 : {nodes_vlos2} # Nodes in vlos\n")
        f.write(f"nodes_gamma2 : {nodes_gamma2} # Nodes in gamma\n")
        f.write(f"nodes_phi2 : {nodes_phi2} # Nodes in phi\n")
        f.write(f"mu_cos : {mu_cos} # mu = cos theta\n")
        f.write(f"abundance : {abundance} # Abundance file\n")
        f.write(f"gas_pressure : {gas_pressure} # Gas Pressure Boundary condition\n")
        f.write(f"# \n")
        f.write(f"# Radomisation Settings\n")
        f.write(f"# \n")
        f.write(f"random_guess : {random_guess} # Create random guesses, 0 = use model as guess\n")
        f.write(f"random_pars : {random_pars} # Randomise these parameters in the file(s) below\n")
        f.write(f"lim_B1 : {lim_B1} # Limits 1 for the randomisation in B in G\n")
        f.write(f"lim_vlos1 : {lim_vlos1} # Limits 1 for the randomisation in vlos in cm/s\n")
        f.write(f"lim_gamma1 : {lim_gamma1} # Limits 1 for the randomisation in the inclination in deg\n")
        f.write(f"lim_azimuth1 : {lim_azimuth1} # Limits 1 for the randomisation in the azimuth in deg\n")
        f.write(f"lim_B2 : {lim_B2} # Limits 2 for the randomisation in B in G\n")
        f.write(f"lim_vlos2 : {lim_vlos2} # Limits 2 for the randomisation in vlos in cm/s\n")
        f.write(f"lim_gamma2 : {lim_gamma2} # Limits 2 for the randomisation in the inclination in deg\n")
        f.write(f"lim_azimuth2 : {lim_azimuth2} # Limits 2 for the randomisation in the azimuth in deg")

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	mode = input ("Which mode do you want to use? [1C/2C/MC]: ")
	if mode == "1C":
		config_1C()
	elif mode == "2C":
		config_2C()
	elif mode == "MC":
		config_MC()
	else:
		print("Mode unknown")

