"""
Changing physical parameters linearly in B, inc, vlos and with splines in T multiple times. Note, that at least log tau
1 to -4 must be covered.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
from scipy.interpolate import PchipInterpolator as Pchip
from shutil import which
import sys
import os
from os.path import exists
import sir
import model_1C as m
import definitions as d


def split_to_float(string, letter=","):
	"""
	Splits the string into a list and converts the elements to floats

	Parameter
	---------
	string : str
		String in the format num1,num2,num3,...
	letter : char
		Defines the separation character
	
	Return
	------
	list containing the information from the string as a float
	"""
	strings = string.split(letter)
	return np.array([float(i) for i in strings])


def help():
	"""
	Information on the function
	"""
	print("create_models - Changes components of a model with two or three nodes")
	print("Usage: python create_models [OPTION]")
	print()
	sir.option("[1. Pos]","Config File")
	sir.option("-save","Savepath for the histograms.")

	sys.exit()


def create_models(conf: dict) -> None:
	"""
	Creates random models depending on the configuration
	
	Parameter
	---------
	config : dict
		Dictionary with the info from the config file

	"""
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/mml.mplstyle'):
		plt.style.use(dirname + '/mml.mplstyle')
		# if dvipng is not installed, don't use latex
		if which('dvipng') is None:
			plt.rcParams["text.usetex"] = "False"
			plt.rcParams["font.family"] = 'sans-serif'
			plt.rcParams["mathtext.fontset"] = 'dejavuserif'
	elif "mml" in plt.style.available:
		plt.style.use('mml')
		# if dvipng is not installed, don't use latex
		if which('dvipng') is None:
			plt.rcParams["text.usetex"] = "False"
			plt.rcParams["font.family"] = 'sans-serif'
			plt.rcParams["mathtext.fontset"] = 'dejavuserif'
	else:
		plt.rcParams["savefig.format"] = "pdf"
	###############################
	# Define variables from input #
	###############################
	path = conf["path"]
	Input = os.path.join(path, conf["model"])
	Output = os.path.join(path, conf["model_out"])
	num = conf["num"]  # Number of random models
	model_nodes = int(conf["model_nodes"])  # If cubic or linear is used

	savepath = path + "/"
	if "-save" in sys.argv:
		savepath = os.path.join(path, sys.argv[sys.argv.index("-save") + 1])
		if not exists(savepath):
			os.makedirs(savepath, exist_ok=True)
	model_pars = conf['model_pars']  # Parameters which should be randomized
	bT = bPe = bvmicro = bB = bvlos = binc = bazi = bz = bPg = brho = False
	if 'T' in model_pars:
		bT = True
	if 'Pe' in model_pars:
		bPe = True
	if 'vmicro' in model_pars:
		bvmicro = True
	if 'B' in model_pars:
		bB = True
	if 'vlos' in model_pars:
		bvlos = True
	if 'gamma' in model_pars:
		binc = True
	if 'phi' in model_pars:
		bazi = True
	if 'z' in model_pars:
		bz = True
	if 'Pg' in model_pars:
		bPg = True
	if 'rho' in model_pars:
		brho = True

	#############
	# LOAD DATA #
	#############
	File = np.loadtxt(Input, skiprows=1)

	# Create arrays with all the columns as rows
	File_T = File.transpose()
	log_tau0 = File_T[0]
	T0 = File_T[1]
	Pe0 = File_T[2]
	vmicro0 = File_T[3]
	B0 = File_T[4]
	vlos0 = File_T[5]
	inc0 = File_T[6]
	azimuth0 = File_T[7]
	if len(File_T) > 8:
		z0 = File_T[8]
		Pg0 = File_T[9]
		rho0 = File_T[10]

	##########################################
	# PREPARATION AND CHANGING OF PARAMETERS #
	##########################################
	create_B = np.array([split_to_float(i, letter=",") for i in conf['create_B'].split(';')])
	create_vlos = np.array([split_to_float(i, letter=",") for i in conf['create_vlos'].split(';')])
	create_gamma = np.array([split_to_float(i, letter=",") for i in conf['create_gamma'].split(';')])
	create_phi = np.array([split_to_float(i, letter=",") for i in conf['create_phi'].split(';')])
	create_points = split_to_float(conf['create_points'])

	model = m.Model()
	model.set_dim(int(num), 1, len(log_tau0))
	model.load = True
	#################
	# LINEAR MODELS #
	#################
	if model_nodes == 2:
		print("-------> Create models with 2 nodes")
		if len(create_points) != 2:
			print("[create_models] The parameter 'create_points' in the config does not have exactly two elements!")
		# Perform 'num' times
		B_0 = np.zeros(num)
		inc_0 = np.zeros(num)
		azi_0 = np.zeros(num)
		vlos_0 = np.zeros(num)
		for i in range(num):
			# Set values to initial model
			model.log_tau[i,0] = np.copy(log_tau0)
			model.T[i,0] = np.copy(T0)
			model.Pe[i,0] = np.copy(Pe0)
			model.vmicro[i,0] = np.copy(vmicro0)
			model.B[i,0] = np.copy(B0)
			model.vlos[i,0] = np.copy(vlos0)
			model.gamma[i,0] = np.copy(inc0)
			model.phi[i,0] = np.copy(azimuth0)
			if len(File_T) > 8:
				model.z[i,0] = np.copy(z0)
				model.Pg[i,0] = np.copy(Pg0)
				model.rho[i,0] = np.copy(rho0)

			######################
			# NEW MAGNETIC FIELD #
			######################
			if bB:
				B_0[i] = np.random.uniform(create_B[0][0], create_B[0][1])  # @ first point
				B_m5 = np.random.uniform(create_B[1][0], create_B[1][1])  # @ second points
				model.B[i,0] = B_0[i] + (B_0[i] - B_m5) / (create_points[1] - create_points[0]) * log_tau0
				B00 = B_0[i]
			else:
				B00 = model.B[i,0][0]
			###################
			# NEW TEMPERATURE #
			###################
			if bT:
				#############################################################################################
				# The two models hsra and cool11 (Collados M., Martínez Pillet V., Ruiz Cobo B., 		  	#
				# Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 (Umbral model for a big spot)) are	#
				# considered. The minimum is the cool11 model, and then I add a factor of the HSRA model.	#
				# The structure of this factor is the following:											#
				# The factors are chosen because of the following thoughts/points:							#
				# - The temperature starts at the model cool11.												#
				# - I create a random factor between 0 and 1.												#
				# - 0 = cool11 Model																		#
				# - 1 = HSRA Model																			#
				# - I implemented some dependency on the magnetic field: If the magnetic field is strong, 	#
				#   the range for the factor is smaller														#
				#############################################################################################

				# Values from HSRA and cool11
				log_taus = np.array([1.4, 1.3, 1.2, 1.1, 1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4,
									0.3, 0.2, 0.1, 0., -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7,
									-0.8, -0.9, -1., -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8,
									-1.9, -2., -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9,
									-3., -3.1, -3.2, -3.3, -3.4, -3.5, -3.6, -3.7, -3.8, -3.9, -4.])

				HSRA_T = np.array([9560., 9390., 9220., 9050., 8880., 8710., 8520., 8290.,
									8030., 7750., 7440., 7140.9, 6860., 6610., 6390., 6200.,
									6035., 5890., 5765., 5650., 5540., 5430., 5330., 5240.,
									5160., 5080., 5010., 4950., 4895., 4840., 4790., 4750.,
									4720., 4690., 4660., 4630., 4600., 4575., 4550., 4525.,
									4490., 4460., 4430., 4405., 4380., 4355., 4330., 4305.,
									4280., 4250., 4225., 4205., 4190., 4175., 4170.])
				cool11_T = np.array([6780.3, 6536., 6291.9, 6048.5, 5806.5, 5569.5, 5340.7, 5117.3,
										4902.9, 4700.4, 4513.9, 4342.3, 4188.8, 4053.4, 3940.5, 3854.,
										3785., 3726.8, 3676.7, 3633.6, 3597.9, 3564.7, 3534.9, 3511.6,
										3498., 3489.4, 3482.2, 3475.6, 3468.9, 3461.6, 3453.6, 3445.2,
										3436.4, 3427., 3417.1, 3406.5, 3395.3, 3383.4, 3370.8, 3357.9,
										3345.1, 3332.4, 3319.2, 3305.5, 3291.1, 3276., 3260.1, 3243.5,
										3225.9, 3207.5, 3188.5, 3170.5, 3155.7, 3142.8, 3129.7]
									)
				HSRA_T -= cool11_T  # Make it relative to cool11 so that with a fac of 1, I get the HSRA model

				# Factor for adding HSRA depending on the magnetic field
				if B00 > 5000:
					factor = np.random.uniform(0.0, 0.6)
				elif B00 > 4000:
					factor = np.random.uniform(0.0, 0.7)
				elif B00 > 3000:
					factor = np.random.uniform(0.0, 0.8)
				elif B00 > 2000:
					factor = np.random.uniform(0.0, 0.9)
				else:
					factor = np.random.uniform(0.6, 1.1)

				# Little perturbation for cool11 model
				cool11_T = cool11_T * np.random.uniform(0.95, 1.05)

				# Add the values
				Ts = cool11_T + factor * HSRA_T

				# Add (only in creating models) additional perturbation in a resulting rotation around log tau -1
				factor = np.random.uniform(0.9, 1.1)
				Ts[Ts > -1] = Ts[Ts > -1] * factor
				Ts[Ts <= -1] = Ts[Ts <= -1] / factor

				model.T[i,0] = np.interp(model.log_tau[i,0], np.flip(log_taus), np.flip(Ts))

			#########################
			# NEW ELECTRON PRESSURE #
			#########################
			if bPe:
				print("Pe randomisation is not implemented")

			#######################
			# NEW MICROTURBULENCE #
			#######################
			if bvmicro:
				print("vmicro randomisation is not implemented")

			############
			# NEW VLOS #
			############
			if bvlos:
				vlos_0[i] = np.random.uniform(create_vlos[0][0], create_vlos[0][1])  # @ 1st point
				vlos_5 = np.random.uniform(create_vlos[1][0], create_vlos[1][1])  # @ 2nd point
				model.vlos[i,0] = vlos_0[i] + (vlos_0[i] - vlos_5) / (create_points[1] - create_points[0]) * log_tau0

			###################
			# NEW INCLINATION #
			###################
			if binc:
				inc_0[i] = np.random.uniform(create_gamma[0][0], create_gamma[0][1])  # @ 1st point
				inc_5 = np.random.uniform(create_gamma[1][0], create_gamma[1][1])  # @ 2nd point
				model.gamma[i,0] = inc_0[i] + (inc_0[i] - inc_5) / (create_points[1] - create_points[0]) * log_tau0

			###############
			# NEW AZIMUTH #
			###############
			if bazi:
				azi_0[i] = np.random.uniform(create_phi[0][0], create_phi[0][1])  # @ 1st point
				azi_5 = np.random.uniform(create_phi[1][0], create_phi[1][1])  # @ 2nd point
				model.phi[i,0] = azi_0[i] + (azi_0[i] - azi_5) / (create_points[1] - create_points[0]) * log_tau0

			##############
			# NEW HEIGHT #
			##############
			if bz:
				print("z randomisation is not implemented")

			####################
			# NEW GAS PRESSURE #
			####################
			if bPg:
				print("Pg randomisation is not implemented")

			###############
			# NEW DENSITY #
			###############
			if brho:
				print("rho randomisation is not implemented")

		model.save(os.path.join(path, Output))

		#######################
		# PLOTTING HISTOGRAMS #
		#######################
		if bB:
			fig, ax = plt.subplots()
			plt.title(r"Histogram of randomly generated magnetic fields @ $\log \tau = $" + str(d.create_points2[0]),
						fontsize=20)
			ax.hist(B_0, bins=20)
			ax.set_xlabel("B [G]")
			ax.set_ylabel("Entries")

			plt.savefig(savepath + "hist_B")

		if bvlos:
			fig, ax = plt.subplots()
			plt.title(
				r"Histogram of randomly generated line of sight velocities @ $\log \tau = $" + str(d.create_points2[0]),
				fontsize=20)
			ax.hist(vlos_0 / 1e5, bins=20)
			ax.set_xlabel(r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}} \right]$")
			ax.set_ylabel("Entries")

			plt.savefig(savepath + "hist_vlos")

		if binc:
			fig, ax = plt.subplots()
			plt.title(r"Histogram of randomly generated inclinations @ $\log \tau = $" + str(d.create_points2[0]),
					fontsize=20)
			ax.hist(inc_0, bins=20)
			ax.set_xlabel(r"$\gamma$ [deg]")
			ax.set_ylabel("Entries")

			plt.savefig(savepath + "hist_inc")

		if bazi:
			fig, ax = plt.subplots()
			plt.title(r"Histogram of randomly generated azimuths @ $\log \tau = $" + str(d.create_points2[0]),
					fontsize=20)
			ax.hist(azi_0, bins=20)
			ax.set_xlabel(r"$\phi$ [deg]")
			ax.set_ylabel("Entries")

			plt.savefig(savepath + "hist_azi")

	elif model_nodes == 3:
		print("-------> Create models with three nodes")
		if len(create_points) != 3:
			print("[create_models] The parameter 'create_points' in the config does not have exactly three elements!")
		B_p1 = np.zeros(num)
		B_m1 = np.zeros(num)
		B_m5 = np.zeros(num)
		inc_p1 = np.zeros(num)
		inc_m1 = np.zeros(num)
		inc_m5 = np.zeros(num)
		azi_p1 = np.zeros(num)
		azi_m1 = np.zeros(num)
		azi_m5 = np.zeros(num)
		vlos_p1 = np.zeros(num)
		vlos_m1 = np.zeros(num)
		vlos_m5 = np.zeros(num)
		# Perform 'num' times
		for i in range(num):
			# Set values to initial model
			model.log_tau[i,0] = np.copy(log_tau0)
			model.T[i,0] = np.copy(T0)
			model.Pe[i,0] = np.copy(Pe0)
			model.vmicro[i,0] = np.copy(vmicro0)
			model.B[i,0] = np.copy(B0)
			model.vlos[i,0] = np.copy(vlos0)
			model.gamma[i,0] = np.copy(inc0)
			model.phi[i,0] = np.copy(azimuth0)
			if len(File_T) > 8:
				model.z[i,0] = np.copy(z0)
				model.Pg[i,0] = np.copy(Pg0)
				model.rho[i,0] = np.copy(rho0)

			######################
			# NEW MAGNETIC FIELD #
			######################
			B00 = B0[0]
			if bB:
				B_p1[i] = np.random.uniform(create_B[0][0], create_B[0][1])  # 1st point
				B_m1[i] = np.random.uniform(B_p1[i] * 0.5, B_p1[i] * 1.1)  # 2nd point

				# For strong magnetic fields => Limit the decreasing part
				if B_p1[i] > 3000:
					# better to not go to too low values to limit the decreasing
					if create_B[1][0] < 500:
						create_B[1][0] = 500
					if create_B[1][1] < 500:
						create_B[1][1] = 1000
				elif B_p1[i] > 2000:
					# better to not go to too low values to limit the decreasing
					if create_B[1][0] < 250:
						create_B[1][0] = 250
					if create_B[1][1] < 600:
						create_B[1][1] = 600

				B_m5[i] = np.random.uniform(create_B[1][0], create_B[1][1])  # @ 3rd point

				spline = inter.CubicSpline(create_points, [B_m5[i], B_m1[i], B_p1[i]], bc_type='natural')
				model.B[i,0] = spline(model.log_tau[i,0])
				B00 = B_p1[i]

			###################
			# NEW TEMPERATURE #
			###################
			if bT:
				###########################################################################################
				# The two models hsra and cool11 (Collados M., Martínez Pillet V., Ruiz Cobo B., 		#
				# Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 (Umbral model for a big spot)) are	#
				# considered. The minimum is the cool11 model, and then I add a factor of the HSRA model.	#
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
				log_taus = np.array([1.4, 1.3, 1.2, 1.1, 1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4,
									0.3, 0.2, 0.1, 0., -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7,
									-0.8, -0.9, -1., -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8,
									-1.9, -2., -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9,
									-3., -3.1, -3.2, -3.3, -3.4, -3.5, -3.6, -3.7, -3.8, -3.9, -4.])

				HSRA_T = np.array([9560., 9390., 9220., 9050., 8880., 8710., 8520., 8290.,
									8030., 7750., 7440., 7140.9, 6860., 6610., 6390., 6200.,
									6035., 5890., 5765., 5650., 5540., 5430., 5330., 5240.,
									5160., 5080., 5010., 4950., 4895., 4840., 4790., 4750.,
									4720., 4690., 4660., 4630., 4600., 4575., 4550., 4525.,
									4490., 4460., 4430., 4405., 4380., 4355., 4330., 4305.,
									4280., 4250., 4225., 4205., 4190., 4175., 4170.]
				)
				cool11_T = np.array([6780.3, 6536., 6291.9, 6048.5, 5806.5, 5569.5, 5340.7, 5117.3,
										4902.9, 4700.4, 4513.9, 4342.3, 4188.8, 4053.4, 3940.5, 3854.,
										3785., 3726.8, 3676.7, 3633.6, 3597.9, 3564.7, 3534.9, 3511.6,
										3498., 3489.4, 3482.2, 3475.6, 3468.9, 3461.6, 3453.6, 3445.2,
										3436.4, 3427., 3417.1, 3406.5, 3395.3, 3383.4, 3370.8, 3357.9,
										3345.1, 3332.4, 3319.2, 3305.5, 3291.1, 3276., 3260.1, 3243.5,
										3225.9, 3207.5, 3188.5, 3170.5, 3155.7, 3142.8, 3129.7]
				)
				HSRA_T -= cool11_T  # Make it relative to cool11 so that with a fac of 1, I get the HSRA model

				# Factor for adding hsra depending on the magnetic field
				if B00 > 5000:
					factor = np.random.uniform(0.0, 0.5)
				elif B00 > 4000:
					factor = np.random.uniform(0.0, 0.6)
				elif B00 > 3000:
					factor = np.random.uniform(0.0, 0.7)
				elif B00 > 2000:
					factor = np.random.uniform(0.0, 0.8)
				else:
					factor = np.random.uniform(0.5, 1.1)

				# Little perturbation for cool11 model
				cool11_T = cool11_T * np.random.uniform(0.95, 1.05)

				# Add the two models
				Ts = cool11_T + factor * HSRA_T

				# Add (only in creating models) additional perturbation in a resulting rotation around log tau -1
				factor = np.random.uniform(0.95, 1.1)

				Ts = Ts * np.linspace(factor, 1 / factor, len(log_taus))

				model.T[i,0] = np.interp(model.log_tau[i,0], np.flip(log_taus), np.flip(Ts))

			#########################
			# NEW ELECTRON PRESSURE	#
			#########################
			if bPe:
				print("Pe randomisation is not implemented")

			#######################
			# NEW MICROTURBULENCE #
			#######################
			if bvmicro:
				print("vmicro randomisation is not implemented")

			############
			# NEW VLOS #
			############
			if bvlos:
				vlos_p1[i] = np.random.uniform(create_vlos[0][0], create_vlos[0][1])  # 1st point
				vlos_m5[i] = np.random.uniform(create_vlos[1][0], create_vlos[1][1])  # 3rd point
				vlos_m1[i] = np.random.uniform(vlos_p1[i], vlos_m5[i])  # 2nd point

				spline = inter.CubicSpline(create_points, [vlos_m5[i], vlos_m1[i], vlos_p1[i]], bc_type='natural')
				model.vlos[i,0] = spline(model.log_tau[i,0]) / 1e5

			###################
			# NEW INCLINATION #
			###################
			if binc:
				inc_p1[i] = np.random.uniform(create_gamma[0][0], create_gamma[0][1])  # 1st point
				inc_m5[i] = np.random.uniform(create_gamma[1][0], create_gamma[1][1])  # 3rd
				inc_m1[i] = np.random.uniform(inc_p1[i], inc_m5[i])  # 2nd point

				# Pchip to prevent overshooting to negative values
				spline = Pchip(create_points, [inc_m5[i], inc_m1[i], inc_p1[i]])
				model.gamma[i,0] = spline(model.log_tau[i,0])

			###############
			# NEW AZIMUTH #
			###############
			if bazi:
				azi_p1[i] = np.random.uniform(create_phi[0][0], create_phi[0][1])  # 1st point
				azi_m5[i] = np.random.uniform(create_phi[1][0], create_phi[1][1])  # 3rd point
				azi_m1[i] = np.random.uniform(azi_p1[i], azi_m5[i])  # 2nd point

				# Pchip to prevent overshooting to negative values
				spline = Pchip(create_points, [azi_m5[i], azi_m1[i], azi_p1[i]])
				model.phi[i,0] = spline(model.log_tau[i,0])

			##############
			# NEW HEIGHT #
			##############
			if bz:
				print("z randomisation is not implemented")

			####################
			# NEW GAS PRESSURE #
			####################
			if bPg:
				print("Pg randomisation is not implemented")

			###############
			# NEW DENSITY #
			###############
			if brho:
				print("rho randomisation is not implemented")

		model.save(os.path.join(path, Output))

		#######################
		# PLOTTING HISTOGRAMS #
		#######################
		if bB:
			fig, ax = plt.subplots()
			plt.title(r"Generated Magnetic Fields", fontsize=20)

			ax.hist(B_p1, bins=20, histtype='step', label=r"@ $\log \tau =  $" + str(create_points[2]))
			ax.hist(B_m1, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[1]))
			ax.hist(B_m5, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[0]))

			ax.set_xlabel("B [G]")
			ax.set_ylabel("Entries")
			ax.legend(loc='upper right')
			plt.savefig(savepath + "hist_B")

		if bvlos:
			fig, ax = plt.subplots()
			plt.title(r"Generated Line of Sight Velocities", fontsize=20)
			ax.hist(vlos_p1 / 1e5, bins=20, histtype='step', label=r"@ $\log \tau =  $" + str(create_points[2]))
			ax.hist(vlos_m1 / 1e5, bins=20, histtype='step', label=r"@ $\log \tau =  $" + str(create_points[1]))
			ax.hist(vlos_m5 / 1e5, bins=20, histtype='step', label=r"@ $\log \tau =  $" + str(create_points[0]))
			ax.set_xlabel(r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}} \right]$")
			ax.set_ylabel("Entries")
			ax.legend(loc='upper right')

			plt.savefig(savepath + "hist_vlos")

		if binc:
			fig, ax = plt.subplots()
			plt.title(r"Generated Inclinations", fontsize=20)
			ax.hist(inc_p1, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[2]))
			ax.hist(inc_m1, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[1]))
			ax.hist(inc_m5, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[0]))

			ax.set_xlabel(r"$\gamma$ [deg]")
			ax.set_ylabel("Entries")
			ax.legend(loc='upper right')

			plt.savefig(savepath + "hist_inc")

		if bazi:
			fig, ax = plt.subplots()
			plt.title(r"Generated Azimuths", fontsize=20)
			ax.hist(azi_p1, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[2]))
			ax.hist(azi_m1, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[1]))
			ax.hist(azi_m5, bins=20, histtype='step', label=r"@ $\log \tau = $" + str(create_points[0]))

			ax.set_xlabel(r"$\phi$ [deg]")
			ax.set_ylabel("Entries")
			ax.legend(loc='upper right')

			plt.savefig(savepath + "hist_azi")

	return


if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	if len(sys.argv) < 2:
		print("[ERROR] No config file provided!")
		sys.exit()
	elif not os.path.exists(sys.argv[1]):
		print("[ERROR] Config file does not exist!")
		sys.exit()
		
	conf = sir.read_config(sys.argv[1])

	create_models(conf)
