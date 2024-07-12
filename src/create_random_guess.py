"""
Create random generated guess models depenidng on the configuration file

Options are 'T', 'Pe', 'vmicro', 'B', 'vlos', 'gamma', 'phi', 'z', 'Pg', 'rho'
"""

import numpy as np 
import sys 
import sir
import definitions as d
import model_atm as m

def help():
	print("create_random_guess - Creates random guesses")
	print("Usage: python create_random_guess.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Config File")
	sir.option("[2. Pos.]","Output")
	sys.exit(1)


def __split_to_float(string, letter = ","):
	"""
	Splits the string into a list and converts the elements to floats

	Parameter
	---------
	string : string
		String in the format num1,num2,num3,...
	letter : char
		Defines the separation character
	
	Return
	------
	list containing the information from the string as a float
	"""

	strings = string.split(letter)
	return np.array([float(i) for i in strings])


def create_guess(model, random_pars, lim_B, lim_vlos, lim_gamma, lim_phi):
	"""
	Create random guesses.

	Parameters
	----------
	model : str
		File name of the base model
	random_pars : list
		List with what is randomised
	lim_B : list
		List with two floats defining the magnetic field at two log tau values
	lim_vlos : list
		List with two floats defining the line-of-sight-velocity at two log tau values
	lim_gamma : list
		List with two floats defining the inclination at two log tau values
	lim_phi : list
		List with two floats defining the azimuth at two log tau values
		
	Return
	------
	Class Model with the random Model
	"""
	
	###################################################################################
	#					Define variables from input					    #
	###################################################################################
	cT = cPe = cvmicro = cB = cvlos = cinc = cazi = cz = cPg = crho = False
	if 'T' in random_pars:
		cT = True
	if 'Pe' in random_pars:
		cPe = True
	if 'vmicro' in random_pars:
		cvmicro = True
	
	if 'B' in random_pars:
		cB = True
		if lim_B is None:
			print("[create_guess] Did you define lim_B in the config file?")
			sys.exit()
	
	if 'vlos' in random_pars:
		cvlos = True
		if lim_vlos is None:
			print("[create_guess] Did you define lim_vlos in the config file?")
			sys.exit()
	
	if 'gamma' in random_pars:
		cinc = True
		if lim_gamma is None:
			print("[create_guess] Did you define lim_gamma in the config file?")
			sys.exit()
	
	if 'phi' in random_pars:
		cazi = True
		if lim_phi is None:
			print("[create_guess] Did you define lim_phi in the config file?")
			sys.exit()

	if 'z' in random_pars:
		cz = True
	if 'Pg' in random_pars:
		cPg = True
	if 'rho' in random_pars:
		crho = True

	#############
	# LOAD DATA #
	#############
	File   = np.loadtxt(model, skiprows=1)
	header = np.genfromtxt(model,max_rows=1)
	File_T = File.transpose()

	mod = m.model_atm(nx = 1, ny = 1, nval=len(File_T[0]))

	#############################
	# SAVE DATA FROM BASE MODEL #
 	#############################
	mod.tau	= File_T[0]
	mod.T[0,0] 	      	= File_T[1]
	mod.Pe[0,0] 	    = File_T[2]
	mod.vmicro[0,0]	= File_T[3]
	mod.B[0,0]			= File_T[4]
	mod.vlos[0,0]		= File_T[5] / 1e5
	mod.gamma[0,0]		= File_T[6]
	mod.phi[0,0]		= File_T[7]
	if len(File_T) > 8:
		mod.z[0,0]		= File_T[8]
		mod.Pg[0,0]		= File_T[9]
		mod.rho[0,0]	= File_T[10]
		mod.full = True

	# Assign Header Information
	mod.vmacro[0,0] = header[0]
	mod.fill[0,0] = header[1]
	mod.stray_light[0,0] = header[2]

	mod.load = True

	###############################
	#	    NEW MAGNETIC FIELD	  #
	###############################
	if cB:
		B0   = np.random.uniform(lim_B[0],lim_B[1]) # @ log tau 1
		mod.B[0,0] = np.ones(shape=mod.B[0,0].shape)*B0 #spline(log_tau)
	else:
		B0 = mod.B[0,0,0]
	#############################
	#	    NEW TEMPERATURE		#
	#############################
	if cT:
		###########################################################################################
		# The two models hsra and cool11 (Collados M., Martínez Pillet V., Ruiz Cobo B., 		  #
		# Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 (Umbral model for a big spot) are  #
		# considered. The minimum is the cool11 model and then I add a factor of the HSRA model.  #
		# The structure of this factor is the following:									      #
		# The factors are chosen because of the following thoughts/points:					      #
		# - The temperature starts at the model cool11.									          #
		# - I create a random factor between 0 and 1.									          #
		# - 0 = cool Model																		  #
		# - 1 = HSRA Model														  	              #
		# - I implemented some dependency on the magnetic field: If the magnetic field is strong, #
		#   the range for the factor is smaller											          #
		###########################################################################################
		# Values from HSRA and cool
		log_taus = d.log_taus
		HSRA_T = np.copy(d.upper_T)
		cool_T = np.copy(d.lower_T)
		
		HSRA_T -= cool_T # Make it relative to cool11 so that with a fac of 1, I get the HSRA model

		# Factor for adding hsra depending on the magnetic field
		factor = np.random.uniform(d.lower_T,d.upper_T)
		# Apply restrictions for stronger magnetic fields
		for i in range(len(d.temp_B)):
			if B0 > d.temp_B[i]:
				factor = np.random.uniform(d.temp_f[i][0], d.temp_f[i][1])
				break
		
		
		# Little perturbation for cool model
		cool_T = cool_T * np.random.uniform(1-d.multiplicative_T, 1+d.multiplicative_T)
		
		Ts = cool_T + factor * HSRA_T
		mod.T[0,0] = np.interp(mod.tau, np.flip(log_taus), np.flip(Ts))
		
		
	#########################
	# NEW ELECTRON PRESSURE #
	#########################
	if cPe:
		print("Pe change not implemented")
	
	#######################
	# NEW MICROTURBULENCE #
	#######################
	if cvmicro:
		print("vmicro change not implemented")

	#############################
	#		    NEW VLOS		#
	#############################
	if cvlos:
		vlos0   = np.random.uniform(lim_vlos[0],lim_vlos[1]) # @ log tau 1
		mod.vlos[0,0] = np.ones(shape=mod.vlos[0,0].shape)*vlos0 / 1e5
	
	#############################
	#	    NEW INCLINATION		#
	#############################
	if cinc:
		inc0   = np.random.uniform(lim_gamma[0],lim_gamma[1]) # @ log tau 1
		mod.gamma[0,0] = np.ones(shape=mod.gamma[0,0].shape)*inc0
	
	#############################
	#	    NEW AZIMUTH		    #
	#############################
	if cazi:
		azi0   = np.random.uniform(lim_phi[0],lim_phi[1])
		mod.phi[0,0] = np.ones(shape=mod.phi[0,0].shape)*azi0
	
	####################################
	#	    NEW HEIGHT			#
	####################################
	if cz:
		print("z change not implemented")
	
	####################################
	#	    NEW GAS PRESSURE	    #
	####################################
	if cPg:
		print("Pg change not implemented")
	
	####################################
	#	    NEW DENSITY		    #
	####################################
	if crho:
		print("rho change not implemented")

	return mod

def create_small_guess(mod, random_pars, factor):
	r"""
	Create random guesses in a small area around the provided values.

	Parameters
	----------
	model : model_atm
		Class with the guesses for the full map
	random_pars : list
		List with what is randomised
	factor : float
		Factor which is used to create a guess around the parameter in the guess file
		E.g. factor = 0.05 would create a parameter space of $[-0.95\cdot P, 1.05\cdot P)$, where P is a model parameter
		
	Return
	------
	Class Model with the random Model

	"""
	
	###################################################################################
	#					Define variables from input					    #
	###################################################################################
	cT = cPe = cvmicro = cB = cvlos = cinc = cazi = cz = cPg = crho = False
	if 'T' in random_pars:
		cT = True
	if 'Pe' in random_pars:
		cPe = True
	if 'vmicro' in random_pars:
		cvmicro = True
	if 'B' in random_pars:
		cB = True
	if 'vlos' in random_pars:
		cvlos = True
	if 'gamma' in random_pars:
		cinc = True
	if 'phi' in random_pars:
		cazi = True
	if 'z' in random_pars:
		cz = True
	if 'Pg' in random_pars:
		cPg = True
	if 'rho' in random_pars:
		crho = True

	for x in range(mod.nx):
		for y in range(mod.ny):		
			###############################
			#	    NEW MAGNETIC FIELD	  #
			###############################
			if cB:
				mod.B[x,y] *= np.random.uniform(1-factor,1+factor)
			
			#############################
			#	    NEW TEMPERATURE		#
			#############################
			if cT:
				mod.T[x,y] *= np.random.uniform(1-factor,1+factor)
				
				
			#########################
			# NEW ELECTRON PRESSURE #
			#########################
			if cPe:
				mod.Pe[x,y] *= np.random.uniform(1-factor,1+factor)
			
			#######################
			# NEW MICROTURBULENCE #
			#######################
			if cvmicro:
				mod.vmicro[x,y] *= np.random.uniform(1-factor,1+factor)

			#############################
			#		    NEW VLOS		#
			#############################
			if cvlos:
				mod.vlos[x,y] *= np.random.uniform(1-factor,1+factor)
			
			#############################
			#	    NEW INCLINATION		#
			#############################
			if cinc:
				mod.gamma[x,y] *= np.random.uniform(1-factor,1+factor)
			
			#############################
			#	    NEW AZIMUTH		    #
			#############################
			if cazi:
				mod.phi[x,y] *= np.random.uniform(1-factor,1+factor)
			
			#############################
			#	    NEW HEIGHT			#
			#############################
			if cz:
				print("z change not implemented")
			
			####################################
			#	    NEW GAS PRESSURE	    #
			####################################
			if cPg:
				print("Pg change not implemented")
			
			####################################
			#	    NEW DENSITY		    #
			####################################
			if crho:
				print("rho change not implemented")

	return mod


def create_guesses_1c(conf, output = "./", number = 0):
	"""
	Create random guess.

	Parameters
	----------
	conf : dict
		Dictionary with the information from the config file
	output : string, optional
		Output path of the models. DefaultL: "./"
	number : int, optional
		added number to the written model. Default: 0

	Return
	------
	None
	"""
	###################################################################################
	#					Define variables from input					    #
	###################################################################################
	if "lim_B" in conf:
		lim_B = __split_to_float(conf["lim_B"])
	else:
		lim_B = None
	if "lim_vlos" in conf:
		lim_vlos = __split_to_float(conf["lim_vlos"])
	else:
		lim_vlos = None
	if "lim_gamma" in conf:
		lim_gamma = __split_to_float(conf["lim_gamma"])
	else:
		lim_gamma = None
	if "lim_phi" in conf:
		lim_phi = __split_to_float(conf[f"lim_phi"])
	else:
		lim_phi = None
	mod = create_guess(conf["model"], conf["random_pars"], lim_B, lim_vlos, lim_gamma, lim_phi)
	# Macroturbulence velocity
	mod.vmacro[0,0] = conf["vmacro"]
	mod.write_model(output + f"{d.model}{number}.mod", 0, 0)


def create_guesses_2c(conf, output = "./", number = 0):
	"""
	Create random guesses.

	Parameters
	----------
	conf : dict
		Dictionary with the information from the config file
	output : string, optional
		Output path of the models. Default: "./"
	number : int, optional
		Added number to the generade output model. Default: 0

	Return
	------
	None
	"""
	###################################################################################
	#					Define variables from input					    #
	###################################################################################

	for j in [1,2]:
		model = conf[f"model{j}"] + str(j)		# Base Model

		if f"lim_B{j}" in conf:
			lim_B = __split_to_float(conf[f"lim_B{j}"])
		else:
			lim_B = None
		if f"lim_vlos{j}" in conf:
			lim_vlos = __split_to_float(conf[f"lim_vlos{j}"])
		else:
			lim_vlos = None
		if f"lim_gamma{j}" in conf:
			lim_gamma = __split_to_float(conf[f"lim_gamma{j}"])
		else:
			lim_gamma = None
		if f"lim_phi{j}" in conf:
			lim_phi = __split_to_float(conf[f"lim_phi{j}"])
		else:
			lim_phi = None

		mod = create_guess(model, conf["random_pars"], lim_B, lim_vlos, lim_gamma, lim_phi)
		mod.fill[0,0] = conf['fill'].split(',')[j-1] # Assign filling factor from config
		
		mod.vmacro[0,0] = conf["vmacro"] # Assign vmacro from config
		if j == 1:
			mod.write_model(output + f"{d.model1}" + str(number) + ".mod", 0, 0)
		else:
			mod.write_model(output + f"{d.model2}" + str(number) + ".mod", 0, 0)

if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])
	
	random_pars = conf['random_pars']
	if "T" in random_pars:
		print("-------> T is changed")
	if "Pe" in random_pars:
		print("-------> Pe is changed")
	if "vmicro" in random_pars:
		print("-------> vmicro is changed")
	if "B" in random_pars:
		print("-------> B is changed")
	if "vlos" in random_pars:
		print("-------> vlos is changed")
	if "gamma" in random_pars:
		print("-------> gamma is changed")
	if "phi" in random_pars:
		print("-------> phi is changed")
	if "z" in random_pars:
		print("-------> z is changed")
	if "Pg" in random_pars:
		print("-------> Pg is changed")
	if "rho" in random_pars:
		print("-------> rho is changed")
	
	if conf["mode"] == "MC" or conf["mode"] == "1C":
		create_guesses_1c(conf, sys.argv[2])
	else:
		create_guesses_2c(conf, sys.argv[2])






	
