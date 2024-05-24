"""

Simulation
==========

This module provides function to create Stokes Profiles by starting with the models. A typical procedure is the following:

 1. Create Models with 1, 2 or 3 Nodes
 2. Perform Synthesis of the Models
 3. Add Noise to the synthesised Stokes Profiles

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
import model_atm as m
import definitions as d
import shutil
import glob
import profile_stk as p




def __split_to_float(string, letter=","):
	"""
	Splits the string into a list and converts the elements to floats

	Parameters
	----------
	string : str
		String in the format num1,num2,num3,...
	letter : char
		Defines the separation character
	
	Returns
	-------
	out : list
		list containing the information from the string as a float
	"""
	strings = string.split(letter)
	return np.array([float(i) for i in strings])

"""
*****************************************************************************
*								ADDING NOISE TO								*
*									PROFILES								*
*****************************************************************************
"""

def add_noise(conf: dict, verbose: bool = True) -> None:
	"""
	Adds noise to the Stokes Profiles

	Parameters
	----------
	conf : dict
		Information from the config file
	verbose : bool, optional
		Print out status. Default: True

	Returns
	-------
	None


	"""
	if verbose:
		print("[STATUS] Add noise to the synthesised profiles")
	##################
	# READ ARGUMENTS #
	##################
	path = conf["path"]
	Input = os.path.join(path, conf["syn_out"] + d.end_stokes)			# Input synthesised model profiles/syn
	Output = os.path.join(path, conf["noise_out"] + d.end_stokes)		# Generic output profiles/noise

	noise_I = conf['noise_I']  # Noise in I
	noise_Q = conf['noise_Q']  # Noise in Q
	noise_U = conf['noise_U']  # Noise in U
	noise_V = conf['noise_V']  # Noise in V

	#############################
	# ADD NOISE TO EACH PROFILE #
	#############################
	syn = p.read_profile(Input)

	noise_Is = np.random.normal(scale=float(noise_I), size=(syn.nx, syn.ny, syn.nw))
	noise_Qs = np.random.normal(scale=float(noise_Q), size=(syn.nx, syn.ny, syn.nw))
	noise_Us = np.random.normal(scale=float(noise_U), size=(syn.nx, syn.ny, syn.nw))
	noise_Vs = np.random.normal(scale=float(noise_V), size=(syn.nx, syn.ny, syn.nw))

	syn.stki = syn.stki + noise_Is
	syn.stkq = syn.stkq + noise_Qs
	syn.stku = syn.stku + noise_Us
	syn.stkv = syn.stkv + noise_Vs

	syn.write(Output)

	return


"""
*****************************************************************************
*								CREATE MODELS								*
*								FOR SYNTHESIS								*
*****************************************************************************
"""
def create_models(conf: dict) -> None:
	"""
	Creates random models depending on the configuration
	
	Parameters
	----------
	config : dict
		Dictionary with the info from the config file

	Returns
	-------
	None

	"""
	dirname = os.path.split(os.path.abspath(__file__))[0]
	plt.rcParams["savefig.format"] = "pdf"
	if d.plt_lib != "":
		plt.style.use(d.plt_lib)
	else:
		if exists(dirname + '/plots/mml.mplstyle'):
			plt.style.use(dirname + '/plots/mml.mplstyle')
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
	###############################
	# Define variables from input #
	###############################
	path = conf["path"]
	Input = os.path.join(path, conf["model"])
	Output = os.path.join(path, conf["model_out"] + d.end_models)
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
	header = np.genfromtxt(Input,max_rows=1)
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
	create_B = np.array([__split_to_float(i, letter=",") for i in conf['create_B'].split(';')])
	create_vlos = np.array([__split_to_float(i, letter=",") for i in conf['create_vlos'].split(';')])
	create_gamma = np.array([__split_to_float(i, letter=",") for i in conf['create_gamma'].split(';')])
	create_phi = np.array([__split_to_float(i, letter=",") for i in conf['create_phi'].split(';')])
	create_points = __split_to_float(conf['create_points'])

	model = m.Model(int(num), 1, len(log_tau0))
	for i in range(num):
		model.vmacro[i,0] = float(conf["vmacro"])
		model.fill[i,0] = header[1]
		model.stray_light[i,0] = header[2]
	model.load = True

	#################
	# Constant MODELS #
	#################
	if model_nodes == 1:
		print("-------> Create models with 1 node")
		if len(create_points) != 2:
			print("[create_models] The parameter 'create_points' in the config does not have exactly two elements!")
		# Perform 'num' times
		B_0 = np.zeros(num)
		inc_0 = np.zeros(num)
		azi_0 = np.zeros(num)
		vlos_0 = np.zeros(num)
		for i in range(num):
			# Set values to initial model
			model.tau = np.copy(log_tau0)
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
				model.B[i,0] = B_0[i] * np.ones(log_tau0.shape)
				B00 = B_0[i]
			else:
				B00 = model.B[i,0][0]
			###################
			# NEW TEMPERATURE #
			###################
			if bT:
				model.T[i,0] = create_temperature(model.tau, B00)

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
				model.vlos[i,0] = vlos_0[i] * np.ones(log_tau0.shape)

			###################
			# NEW INCLINATION #
			###################
			if binc:
				inc_0[i] = np.random.uniform(create_gamma[0][0], create_gamma[0][1])  # @ 1st point
				model.gamma[i,0] = inc_0[i] * np.ones(log_tau0.shape)

			###############
			# NEW AZIMUTH #
			###############
			if bazi:
				azi_0[i] = np.random.uniform(create_phi[0][0], create_phi[0][1])  # @ 1st point
				model.phi[i,0] = azi_0[i] * np.ones(log_tau0.shape)

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

	#################
	# LINEAR MODELS #
	#################
	elif model_nodes == 2:
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
			model.tau = np.copy(log_tau0)
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
				model.T[i,0] = create_temperature(model.tau, B00)

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
			model.tau = np.copy(log_tau0)
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
				model.B[i,0] = spline(model.tau)
				B00 = B_p1[i]

			###################
			# NEW TEMPERATURE #
			###################
			if bT:
				model.T[i,0] = create_temperature(model.tau, B00)

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
				model.vlos[i,0] = spline(model.tau) / 1e5

			###################
			# NEW INCLINATION #
			###################
			if binc:
				inc_p1[i] = np.random.uniform(create_gamma[0][0], create_gamma[0][1])  # 1st point
				inc_m5[i] = np.random.uniform(create_gamma[1][0], create_gamma[1][1])  # 3rd
				inc_m1[i] = np.random.uniform(inc_p1[i], inc_m5[i])  # 2nd point

				# Pchip to prevent overshooting to negative values
				spline = Pchip(create_points, [inc_m5[i], inc_m1[i], inc_p1[i]])
				model.gamma[i,0] = spline(model.tau)

			###############
			# NEW AZIMUTH #
			###############
			if bazi:
				azi_p1[i] = np.random.uniform(create_phi[0][0], create_phi[0][1])  # 1st point
				azi_m5[i] = np.random.uniform(create_phi[1][0], create_phi[1][1])  # 3rd point
				azi_m1[i] = np.random.uniform(azi_p1[i], azi_m5[i])  # 2nd point

				# Pchip to prevent overshooting to negative values
				spline = Pchip(create_points, [azi_m5[i], azi_m1[i], azi_p1[i]])
				model.phi[i,0] = spline(model.tau)

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

		model.write(os.path.join(path, Output))

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
	else:
		print(f"[ERROR] Number of nodes {model_nodes} not implemented!")

	return

def create_temperature(tau, B = 0):
	r"""
	Creates a random temperature in the parameter space.

	Parameters
	----------
	tau : numpy array
		Array with the $\log\tau$ of the model
	B : float, optional
		Value of the magnetic field in the deepest layer. Default: 0

	Returns
	-------
	out : numpy array
		Interpolated numy array with a random temperature atmosphere
	"""
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
	if B > 5000:
		factor = np.random.uniform(0.0, 0.6)
	elif B > 4000:
		factor = np.random.uniform(0.0, 0.7)
	elif B > 3000:
		factor = np.random.uniform(0.0, 0.8)
	elif B > 2000:
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

	return np.interp(tau, np.flip(log_taus), np.flip(Ts))

def synthesis(conf, comm, rank, size, MPI):
	"""
	Performs the synthesis of all the models.

	Parameters
	----------
	config : dict
		Dcitionary with config information
	rank : int
		Number of this process
	comm : Intracomm
		Communicator from MPI
	size : int
		Number of processes
	MPI : library
		Library MPI

	Returns
	-------
	None

	"""
	
	comm.barrier()

	###################################################
	#		READ PARAMETERS FROM CONFIG FILE		#	
	###################################################
	# Define parameters for easier understanding
	path = conf["path"]
	abundance_file = conf['abundance'] # Abundance file	
	
	models = m.read_model(os.path.join(path, conf['model_out'] + d.end_models))

	####################################
	#	CREATE GRID AND CONFIG FILE	#
	####################################
	if rank == 0:
		sir.write_control(os.path.join(path,d.syn_trol_file), conf, Type="syn")
		sir.write_grid(conf,os.path.join(path,d.Grid))
			
	if rank == 0:
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)
	
	comm.barrier() # Continues when all processes entered this function

	
	####################################
	#		START INVERSION PART	#
	####################################

	performed_models = 0 # Counts how many models are performed
	finished_jobs = 0

	tasks = sir.create_task_folder_list(conf['num'])

	for i in range(rank, conf['num'], size):
		# Create task folder
		task_folder = os.path.join(path, tasks['folders'][i])
		os.makedirs(task_folder, exist_ok=True)
		
		# Copy SIR files and guess model
		sir_files = [d.syn_trol_file, "sir.x", conf['line'], d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(task_folder, sir_file))
		
		# Extract model
		models.write_model(os.path.join(task_folder,d.model_syn), i, 0)
	
		# Perform synthesis
		os.chdir(task_folder)
		os.system("echo " + d.syn_trol_file + " | " + " ./sir.x >/dev/null 2>/dev/null")
		
		performed_models += 1

		if finished_jobs < (conf['num'] - conf['num'] % size):
			finished_jobs = comm.allreduce(performed_models, op=MPI.SUM)

		if rank == 0:
			print(f"\rTotal finished Jobs: {finished_jobs}", end='', flush=False)

		os.chdir('../') # Go back in case relative paths are used
	os.chdir(path)
	# Collect data and save it
	if rank == 0:
		print(f"\rTotal finished Jobs: {conf['num']}", end='', flush=False)
		output =  os.path.join(path,conf["syn_out"]) # Output file
		atoms = [i.split(",") for i in conf['atoms']]

		# Read the profiles
		stk = p.Profile(conf['num'],1,0)
		stk.read_results_MC(path, tasks, d.profile)
		
		stk.write(f"{os.path.join(conf['path'],conf['syn_out'] + d.end_models)}")
		
		for i in range(conf['num']):
			shutil.rmtree(tasks['folders'][i])
		print(f"\r-------> Finished with {conf['num']} synthesised models.")

	comm.barrier()
