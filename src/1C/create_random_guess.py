"""
Create random generated guess models depenidng on the configuration file

Options are 'T', 'Pe', 'vmicro', 'B', 'vlos', 'gamma', 'phi', 'z', 'Pg', 'rho'
"""

import numpy as np 
import sys 
import sir
from os.path import exists
import os
import scipy.interpolate as inter
import definitions as d
import matplotlib.pyplot as plt
import model_1C as m

def help():
	print("create_random_guesses - Creates random guesses")
	print("Usage: python create_random_guesses.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Config File")
	sir.option("[2. Pos.]","Output")
	sys.exit(1)


def split_to_float(string, letter = ","):
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


def create_guesses(conf, output = "./", number = -1):
	"""
	Create random guesses.

	Parameters
	----------
	conf : dict
		Dictionary with the information from the config file
	output : string
		Output path of the models
	number : int, optional
		only generate one random guess and add this number to the file name. Default: -1 (not used)

	Return
	------
	None
	"""
	###################################################################################
	#					Define variables from input					    #
	###################################################################################
	model = conf["model"]		# Base Model
	num   = conf['random_guess']  # Number of random guesses
	
	random_pars = conf['random_pars'] # Parameters which should be randomized
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

	lim_B		= split_to_float(conf["lim_B"])
	lim_vlos	= split_to_float(conf["lim_vlos"])
	lim_gamma	= split_to_float(conf["lim_gamma"])
	lim_phi		= split_to_float(conf["lim_phi"])

	###################################################################################
	#						    LOAD DATA							   #
	###################################################################################
	File   = np.loadtxt(model, skiprows=1)
	File_T    = File.transpose()
	Header = np.loadtxt(model, max_rows=1)

	# Create header
	Header = "   " + str(conf["vmacro"]) + "  " + str(1.0) + "  " + str(Header[2])

	###################################################################################
	#			    PREPARATION AND CHANGING OF PARAMETERS				    #
	###################################################################################
	if not exists(output):
		os.mkdir(output)

	# Define the number of guess models
	Numbers = range(num)
	if number != -1:
		Numbers = [number-1]

	mod = m.Model(nx = len(Numbers), ny = 1, npar=len(File_T[0]))

	# Perform 'num' times
	for i in range(len(Numbers)):
		# Arrays with pars from the model
		# Create arrays with all the columns as rows
		
		mod.log_tau[i,0]	= File_T[0]
		mod.T[i,0] 	      	= File_T[1]
		mod.Pe[i,0] 	    = File_T[2]
		mod.vmicro[i,0]	= File_T[3]
		mod.B[i,0]			= File_T[4]
		mod.vlos[i,0]		= File_T[5] / 1e5
		mod.gamma[i,0]		= File_T[6]
		mod.phi[i,0]		= File_T[7]
		if len(File_T) > 8:
			mod.z[i,0]		= File_T[8]
			mod.Pg[i,0]		= File_T[9]
			mod.rho[i,0]	= File_T[10]
		
		####################################
		#	    NEW MAGNETIC FIELD	  #
		####################################
		if cB:
			B0   = np.random.uniform(lim_B[0],lim_B[1]) # @ log tau 1
			mod.B[i,0] = np.ones(shape=mod.B[i,0].shape)*B0 #spline(log_tau)
		else:
			B0 = mod.B[i,0,0]
		####################################
		#	    NEW TEMPERATURE			#
		####################################
		if cT:
			###########################################################################################
			# The two models hsra and cool11 (Collados M., Martínez Pillet V., Ruiz Cobo B., 		#
			# Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 (Umbral model for a big spot) are	#
			# considered. The minimum is the cool11 model and then I add a factor of the HSRA model.	#
			# The structure of this factor is the following:									#
			# The factors are chosen because of the following thoughts/points:					#
			# - The temperature starts at the model cool11.									#
			# - I create a random factor between 0 and 1.									#
			# - 0 = cool11 Model														#
			# - 1 = HSRA Model															#
			# - I implemented some dependency on the magnetic field: If the magnetic field is strong, #
			#   the range for the factor is smaller											#
			###########################################################################################

			# Values from HSRA and cool11
			log_taus = np.array([ 1.4,  1.3,  1.2,  1.1,  1. ,  0.9,  0.8,  0.7,  0.6,  0.5,  0.4,
							   0.3,  0.2,  0.1,  0. , -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7,
							  -0.8, -0.9, -1. , -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8,
							  -1.9, -2. , -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9,
							  -3. , -3.1, -3.2, -3.3, -3.4, -3.5, -3.6, -3.7, -3.8, -3.9, -4. ])

			# 
			HSRA_T = np.array([9560. , 9390. , 9220. , 9050. , 8880. , 8710. , 8520. , 8290. ,
						  8030. , 7750. , 7440. , 7140.9, 6860. , 6610. , 6390. , 6200. ,
						  6035. , 5890. , 5765. , 5650. , 5540. , 5430. , 5330. , 5240. ,
						  5160. , 5080. , 5010. , 4950. , 4895. , 4840. , 4790. , 4750. ,
						  4720. , 4690. , 4660. , 4630. , 4600. , 4575. , 4550. , 4525. ,
						  4490. , 4460. , 4430. , 4405. , 4380. , 4355. , 4330. , 4305. ,
					       4280. , 4250. , 4225. , 4205. , 4190. , 4175. , 4170. ])
			cool_T = np.array([6780.3, 6536. , 6291.9, 6048.5, 5806.5, 5569.5, 5340.7, 5117.3,
						  4902.9, 4700.4, 4513.9, 4342.3, 4188.8, 4053.4, 3940.5, 3854. ,
						  3785. , 3726.8, 3676.7, 3633.6, 3597.9, 3564.7, 3534.9, 3511.6,
						  3498. , 3489.4, 3482.2, 3475.6, 3468.9, 3461.6, 3453.6, 3445.2,
						  3436.4, 3427. , 3417.1, 3406.5, 3395.3, 3383.4, 3370.8, 3357.9,
						  3345.1, 3332.4, 3319.2, 3305.5, 3291.1, 3276. , 3260.1, 3243.5,
						  3225.9, 3207.5, 3188.5, 3170.5, 3155.7, 3142.8, 3129.7]
						)
			HSRA_T -= cool_T # Make it relative to cool11 so that with a fac of 1, I get the HSRA model

			# Factor for adding hsra depending on the magnetic field
			if B0 > 5000:
				factor = np.random.uniform(0.0,0.5)
			elif B0 > 4000:
				factor = np.random.uniform(0.0,0.6)
			elif B0 > 3000:
				factor = np.random.uniform(0.0,0.7)
			elif B0 > 2000:
				factor = np.random.uniform(0.0,0.8)
			else:
				factor = np.random.uniform(0.5,1.1)
			
			# Little perturbation for cool model
			cool_T = cool_T * np.random.uniform(0.95,1.05)
			
			Ts = cool_T + factor * HSRA_T

			mod.T[i,0] = np.interp(mod.log_tau[i,0], np.flip(log_taus), np.flip(Ts))
			
			
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


		#####################################
		#		    NEW VLOS				#
		#####################################
		if cvlos:

			vlos0   = np.random.uniform(lim_vlos[0],lim_vlos[1]) # @ log tau 1
			mod.vlos[i,0] = np.ones(shape=mod.vlos[i,0].shape)*vlos0 / 1e5

		####################################
		#	    NEW INCLINATION			   #
		####################################
		if cinc:
			inc0   = np.random.uniform(lim_gamma[0],lim_gamma[1]) # @ log tau 1
			mod.gamma[i,0] = np.ones(shape=mod.gamma[i,0].shape)*inc0

		####################################
		#	    NEW AZIMUTH		    #
		####################################
		if cazi:
			azi0   = np.random.uniform(lim_phi[0],lim_phi[1])
			mod.phi[i,0] = np.ones(shape=mod.phi[i,0].shape)*azi0
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
	
		# For testing purpose to visualize how it looks (uncomment if wanted)
		"""
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

		ax1.plot(mod.log_tau[i,0],mod.T[i,0])
		ax2.plot(mod.log_tau[i,0],mod.B[i,0])
		ax3.plot(mod.log_tau[i,0],mod.vlos[i,0])
		ax4.plot(mod.log_tau[i,0],mod.gamma[i,0])
		ax1.set_ylabel("T")
		ax2.set_ylabel("B")
		x3.set_ylabel("vlos")
		ax4.set_ylabel(r"$\gamma$")
		ax1.set_xlim(mod.log_tau[i,0,0],mod.log_tau[i,0,-1])
		ax2.set_xlim(mod.log_tau[i,0,0],mod.log_tau[i,0,-1])
		ax3.set_xlim(mod.log_tau[i,0,0],mod.log_tau[i,0,-1])
		ax4.set_xlim(mod.log_tau[i,0,0],mod.log_tau[i,0,-1])

		plt.show()
		"""

	# Write the new models
	for i in range(len(Numbers)):
		mod.write(output + model.replace(".mod","") + str(Numbers[i]+1) + ".mod", Header,i,0)

	

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

	create_guesses(conf, sys.argv[2])






	
