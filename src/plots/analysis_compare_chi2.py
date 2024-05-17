"""
Compares two simulations by computing and plotting the chi2 value as well as computing the mean value.
"""

import numpy as np 
import sys, os
sys.path.append(sys.path[0] + "/..")
import sir
import model_atm as m
import definitions as d
import matplotlib.pyplot as plt
from os.path import exists

def _help():
	"""
	Help Page
	"""
	print("analysis_compare_chi2 - Computes the chi2 of two simulations")
	print("Usage: python analysis_compare_chi2 [OPTION]")
	print()
	sir.option("[1. Pos]","Fitted model 1 (generic, adds a number and .mod)")
	sir.option("[2. Pos]","Actual model 1")
	sir.option("[3. Pos]","Number of models 1")
	sir.option("[4. Pos]","Fitted model 2")
	sir.option("[5. Pos]","Actual model 2")
	sir.option("[6. Pos]","Number of models 2")
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
	sir.option("-T","Plot temperature in K")
	sir.option("-B","Plot magentic field strength in Gauss")
	sir.option("-vlos","Plot line of sight velocity in cm/s")
	sir.option("-inc","Plot inclination by subtracting in deg")
	sir.option("-azi","Plot azimuth by adding in deg")
	sir.option("-add1","Additional text in label for model 1")
	sir.option("-add2","Additional text in label for model 2")
	sys.exit()


def _chi2_err(y_fit, y_obs, yerr):
	r"""
	Computes the merit-function $\chi^2$.
	
	Parameters
	----------
	y_obs : numpy array
		Array with the observed values
	y_fit : numpy array
		Array with computed values, e.g. from a fit
	yerr : numpy.array
		Array containing the errors on the fit/obs

	Return
	------
	out : float
		$\chi^2$-value
	
	"""
	return np.sum((y_fit-y_obs)**2 / yerr**2)

def analysis_compare_chi2(conf1, conf2):
	"""
	Compares the chi2 of two simulations for different physical parameters.
	
	Parameters
	----------
	conf1 : dict
		Configuration for simulation 1
	conf2 : dict
		Configuration for simulation 2

	Returns
	-------
	None

	Other Parameters
	----------------
	Additional parameters given as an argument when the script is executed.
	-save [str], optional
		Additional save path. Default: './'
	-add [str]
		Additional text in filenames
	-label [str]
		Add label text
	-T
		Plot temperature in K
	-B
		Plot magentic field strength in Gauss
	-vlos
		Plot line of sight velocity in cm/s
	-inc
		Plot inclination by subtracting in deg
	-azi
		Plot azimuth by adding in deg
	-add1 [str]
		Additional text in label for model 1
	-add2 [str]
		Additional text in label for model 2


	"""

	# Import library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	plt.rcParams["savefig.format"] = "pdf"
	if d.plt_lib != "":
		plt.style.use(d.plt_lib)
	else:
		if exists(dirname + '/mml.mplstyle'):
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


	###############################################
	# READ INPUT, WRITE DEFINITIONS AND LOAD DATA #
	###############################################
	path1 = conf1["path"]
	num1 = conf1["num"]
	path2 = conf2["path"]
	num2 = conf2["num"]

	fit1 = m.model(os.path.join(path1, conf1["inv_out"]) + d.end_models)
	err1 = m.error(os.path.join(path1, conf1["inv_out"]) + d.end_errors)
	mod1 = m.model(os.path.join(path1, conf1["model_out"]))

	fit2 = m.model(os.path.join(path2, conf2["inv_out"]) + d.end_models)
	err2 = m.error(os.path.join(path2, conf2["inv_out"]) + d.end_errors)
	mod2 = m.model(os.path.join(path2, conf2["model_out"]))

	# Additional savepath
	savepath = ''
	if "-save" in sys.argv:
		savepath = sys.argv[sys.argv.index("-save")+1]

	# Additional text in saved figures
	add_text = ''
	if "-add" in sys.argv:
		add_text = sys.argv[sys.argv.index("-add")+1]


	# Additional labels
	add_label = '_'
	if "-label" in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]
	add1 = '_'
	if "-add1" in sys.argv:
		add1 = " " + sys.argv[sys.argv.index("-add1")+1]
	add2 = '_'
	if len(sys.argv) > 9:
		add2 = " " + sys.argv[sys.argv.index("-add2")+1]

	# Plotting settings
	Markers = ["-", '--', 'dotted', 'dashdotdotted', 'densely dashed']

	linestyle_str = [
		'solid',	 # Same as (0, ()) or '-'
		'dotted',    # Same as (0, (1, 1)) or ':'
		'dashed',    # Same as '--'
		'dashdot',
		'solid',	 # Same as (0, ()) or '-'
		'dotted',    # Same as (0, (1, 1)) or ':'
		'dashed',    # Same as '--'
		'dashdot',
		(0, (3,10,1,10))
		] 

	#######################################
	# DETERMINE WHAT PARAMETER ARE CHOSEN #
	#######################################

	inputs = ["-T", "-B", "-vlos", "-inc", "-azi"]
	att = ["T", "B", "vlos", "gamma", "phi"]
	index = [1, 4, 5, 6, 7]
	ind_used = []  # To determine which values should be computed and plotted
	limits = [2e6, 25000, 8000, 900, 2e5]
	# Labels and titles for the saved plot
	labels = ["T [K]", "B [Gaus]",
			r"$v_{\mathrm{los}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}^{-1}}\right]$",
			r"$\gamma$ [deg]", r"$\phi$ [deg]"]
	titles = [r"Temperature T", r"Magnetic Field Strength B",
			r"Line-of-sight Velocity $v_{\mathrm{los}}$",
			r"Inclination $\gamma$", r"Azimuth $\phi$"]


	# Count how many parameters should be plotted
	counts = 0
	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			ind_used.append(i)
			counts += 1

	##########################################################
	# Compute chi for each model and each physical parameter #
	##########################################################
	# Create the chi 
	chi1 = np.zeros(shape=(counts, num1))
	chi2 = np.zeros(shape=(counts, num2))

	# Go through each model and compute chi2 for each wished physical parameter for model 1
	for i in range(num1):
		# Load fit data
		Tf, Bf, vlosf, incf, azimuthf = fit1.T[:,0,:], fit1.B[:,0,:], fit1.vlos[:,0,:], fit1.gamma[:,0,:], fit1.phi[:,0,:]
		eT, eB, evlos, einc, eazimuth = err1.T[:,0,:], err1.B[:,0,:], err1.vlos[:,0,:], err1.gamma[:,0,:], err1.phi[:,0,:]
		Tm, Bm, vlosm, incm, azimuthm = mod1.T[:,0,:], mod1.B[:,0,:], mod1.vlos[:,0,:], mod1.gamma[:,0,:], mod1.phi[:,0,:]

		# Compute chi2 including the error
		n = 0
		if "-T" in sys.argv:
			chi1[n, i] = _chi2_err(Tf, Tm, eT)
			n += 1
		if "-B" in sys.argv:
			chi1[n, i] = _chi2_err(Bf, Bm, eB)
			n += 1
		if "-vlos" in sys.argv:
			chi1[n, i] = _chi2_err(vlosf, vlosm, evlos)
			n += 1
		if "-inc" in sys.argv:
			chi1[n, i] = _chi2_err(incf, incm, einc)
			n += 1
		if "-azi" in sys.argv:
			chi1[n, i] = _chi2_err(azimuthf, azimuthm, eazimuth)
			n += 1

	# Compute chi2 for second model
	for i in range(num2):
		# Load fit data
		Tf, Bf, vlosf, incf, azimuthf = fit2.T[:,0,:],fit2.B[:,0,:],fit2.vlos[:,0,:],fit2.gamma[:,0,:],fit2.phi[:,0,:]
		eT, eB, evlos, einc, eazimuth = err2.T[:,0,:],err2.B[:,0,:],err2.vlos[:,0,:],err2.gamma[:,0,:],err2.phi[:,0,:]
		Tm, Bm, vlosm, incm, azimuthm = mod2.T[:,0,:],mod2.B[:,0,:],mod2.vlos[:,0,:],mod2.gamma[:,0,:],mod2.phi[:,0,:]

		# Compute chi2 including the error
		n = 0
		if "-T" in sys.argv:
			chi2[n, i] = _chi2_err(Tf, Tm, eT)
			n += 1
		if "-B" in sys.argv:
			chi2[n, i] = _chi2_err(Bf, Bm, eB)
			n += 1
		if "-vlos" in sys.argv:
			chi2[n, i] = _chi2_err(vlosf, vlosm, evlos)
			n += 1
		if "-inc" in sys.argv:
			chi2[n, i] = _chi2_err(incf, incm, einc)
			n += 1
		if "-azi" in sys.argv:
			chi2[n, i] = _chi2_err(azimuthf, azimuthm, eazimuth)
			n += 1

	#######################################################################
	#					   PRINT INFORMATION				    #
	#######################################################################
	# ! Compute the mean and print it out ! #

	mean1 = np.mean(chi1, axis = 1)
	mean2 = np.mean(chi2, axis = 1)

	print()
	print()
	print("		|", end='')
	for i in range(counts):
		size = 5 - len(inputs[ind_used[i]][1:])
		white = "		     "[0:size+2] + "|"
		print("   " , inputs[ind_used[i]][1:], white, end='')

	print()
	print(" Model 1  |", end='')
	for i in range(counts):
		print(" %1.3E  |" % mean1[i], end='')
	print()
	print(" Model 2  |", end='')
	for i in range(counts):
		print(" %1.3E  |" % mean2[i], end='')
	print()
	print()


	if not exists(savepath):
		os.mkdir(savepath)
	# ! Plot physical parameters ! 
	for i in range(counts):
			
			fig, ax1 = plt.subplots(figsize=(8,8))

			plt.title(titles[ind_used[i]])
			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)
			
			# Determine equal bins for both models
			bins = np.linspace(
					np.min(np.append(chi1[i,:],chi2[i,:]))*0.9,
					np.max(np.append(chi1[i,:],chi2[i,:]))*1.1,
					20
				  ) 
			bins = np.linspace(0, limits[i],20)
			ax1.hist(chi1[i,:], label="Model 1" + add1, bins=bins, histtype='step')
			ax1.hist(chi2[i,:], label="Model 2" + add2, bins=bins, histtype='step')

			# Set labels
			ax1.set_xlabel(r"$\chi^2$")
			ax1.set_ylabel(r"Entries")
			
			ax1.legend()

			plt.tight_layout()

			plt.savefig(savepath + "sir_chi2_" + inputs[ind_used[i]][1:] + add_text + ".png")
			

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
	conf1 = sir.read_config(sys.argv[1])
	conf2 = sir.read_config(sys.argv[2])
	if conf1["mode"] == "MC" and conf2['mode']:
		analysis_compare_chi2(conf1,conf2)
	else:
		if conf1["mode"] == "MC":
			print(f"[visualizer] Mode '{conf2['mode']}' of config 2 unknown or not defined")
		elif conf2["mode"] == "MC":
			print(f"[visualizer] Mode '{conf1['mode']}' of config 1 unknown or not defined")
		else:
			print(f"[visualizer] Mode '{conf1['mode']}' and '{conf2['mode']}'unknown or not defined")
		
