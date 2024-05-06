"""

Analyses the quality of the Monte Carlo Simulation by computing the standard deviations

"""

import numpy as np 
import sys, os
sys.path.append(sys.path[0] + "/../..")
sys.path.append(os.path.join(sys.path[0], "/../../tools"))
import sir
import definitions as d
from model import * # Class Model

import matplotlib.pyplot as plt
from os.path import exists
from tools.change_config_path import change_config_path

def _std_dev(data, syn):
	"""
	Computes the standard deviation along log tau

	"""

	std = np.zeros(data.shape[1])
	for j in range(data.shape[1]):
		Sum = 0.0
		for i in range(data.shape[0]):
			Sum += (data[i][j]-syn[i][j])**2
		std[j] = np.sqrt(Sum/(data.shape[0]-1))
	return std


def _help():
	"""
	Help Page

	"""
	print("analysis - Analyses the quality of the Monte Carlo Simulation by plotting the standard deviations")
	print("Usage: python analysis.py [OPTION]")
	print()
	sir.option("[1. Pos]","Config file")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-T","Compare temperature in K")
	sir.option("-B","Compare magentic field strength in Gauss")
	sir.option("-vlos","Compare line of sight velocity in cm/s")
	sir.option("-inc","Compare inclination by subtracting in deg")
	sir.option("-azi","Compare azimuth by adding in deg")
	sir.option("-title","Title of the 4 figures plot")
	sir.option("-xtitle","x position of title in Stokes plot (optional)")
	sir.option("-max","Plot the maximum differences")
	sir.option("-vertical","Plot the last plot vertical")
	sir.option("-v:","Print out tables with values at different log taus.")
	
	sys.exit()

def _get_attribute(model, att):
	'''
	Returns a model parameter depending on a string

	Parameters
	----------
	model : Model
		Class with the Model
	att : str
		String which attribute is returned

	Returns
	-------
	None

	'''
	if att == "tau":
		return model.log_tau
	elif att == "T":
		return model.T[:,0,:]
	elif att == "Pe":
		return model.Pe[:,0,:]
	elif att == "vmicro":
		return model.vmicro[:,0,:]
	elif att == "B":
		return model.B[:,0,:]
	elif att == "vlos":
		return model.vlos[:,0,:]
	elif att == "gamma":
		return model.gamma[:,0,:]
	elif att == "phi":
		return model.phi[:,0,:]
	elif att == "z":
		return model.z[:,0,:]
	elif att == "Pg":
		return model.Pg[:,0,:]
	elif att == "rho":
		return model.rho[:,0,:]
	else:
		return 0

def analysis(conf):
	"""

	Analysis of a MC simulation by plotting the standard deviations

	Parameters
	----------
	config : dict
		Configurations

	Returns
	-------
	None

	Other Parameters
	----------------

	There are optional options which change the plots. Additional parameters given as an argument when the script is executed.

	-save [str], optional
		Additional save path. Default './'.
	-add [str]
		Additional text in filenames.
	-T
		Compare temperature in K
	-B
		Compare magentic field strength in Gauss
	-vlos
		Compare line of sight velocity in cm/s
	-inc
		Compare inclination by subtracting in deg
	-azi
		Compare azimuth by adding in deg
	-title [str]
		Title of the 4 figures plot
	-xtitle [float]
		x position of title in Stokes plot
	-max
		Plot the maximum differences
	-vertical
		Plot the last plot vertical
	-v
		Print out tables with values at different log taus.

	"""

	# Import library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/../../mml.mplstyle'):
		plt.style.use(dirname + '/../../mml.mplstyle')
	elif "mml" in plt.style.available:
		plt.style.use('mml')

	# Check if path exists
	if not exists(conf['path']):
		Inp = input("[NOTE] Path does not exist. You want to overwrite it with the actual path? [y/n] ")
		if Inp == "y":
			change_config_path(conf,os.path.abspath(os.getcwd()))


	#########################################################################
	#	    READ INPUT, WRITE DEFINITIONS AND LOAD DATA				  #
	#########################################################################
	path = conf["path"]
	data = read_model(os.path.join(path,conf["inv_out"]) + d.end_models) # Data from fit
	syn  = read_model(os.path.join(path,conf["model_out"])) # Data

	# Correct phi range
	syn.correct_phi()
	data.correct_phi()

	# Additional labels
	savepath = ''
	if "-save" in sys.argv:
		savepath = path + "/" + sys.argv[sys.argv.index("-save")+1]

	# Additional text in saved figures
	add = ''
	if "-add" in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Title of the last plot
	title = ''
	if "-title" in sys.argv:
		title = sys.argv[sys.argv.index("-title")+1]

	# Title x position
	xtitle = -1
	if '-xtitle' in sys.argv:
		xtitle = float(sys.argv[sys.argv.index("-xtitle")+1])

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
	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]  # Get colors used in the actual cycle

	#######################################################################
	#		    DETERMINE WHAT PARAMETER ARE CHOSEN				#
	#######################################################################

	inputs = ["__","-T", "-Pe", "-vmicro", "-B", "-vlos", "-inc", "-azi", "-z", "-Pg", "-rho"]
	# Labels and titles for the saved plot
	att = ["tau","T", "Pe", "vmicro", "B", "vlos", "gamma", "phi", "z", "Pg", "rho"] # For getting the value from the class
	labels_y = ["",r"$_T$", "", "", r"$_B$", r"$_{\mathrm{v}_{\mathrm{los}}}$", r"$_\gamma$", r"$_\phi$", r"$_z$", r"$_{P_g}$", r"$_{rho}$"]
	labels = ["",r"$_T$ [K]", "", "", r"$_B$ [G]", r"$_{\mathrm{v}_{\mathrm{los}}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$_\gamma$ [deg]", r"$_\phi$ [deg]", r"$_z$ [km]", "$_{P_g}$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$_{\rho}$ $\left[\frac{\mathrm{g}}{\mathrm{cm}^3}\right]$"]

	titles = ["", r"Temperature T", "", "", r"Magnetic Field Strength B",
			  r"Line-of-Sight Velocity $v_{\mathrm{los}}$", r"Inclination $\gamma$",
			  r"Azimuth $\phi$", r"Height z", r"Gas Pressure $P_g$"]


	if savepath != '' and not exists(savepath):
		os.mkdir(savepath)

	#######################################################################
	#					   Plot physical parameters			  #
	#######################################################################
	lim_max = -3
	# Cut data so that the plot limits are adjusted to the shorted range
	syn.set_limit(lim_max)
	data.set_limit(lim_max)
	log_tau = data.log_tau
	lim_max = data.log_tau[-1]

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:

			fig, ax1 = plt.subplots(figsize=(8,7))

			plt.title(titles[i])

			# Standard deviation
			std = _std_dev(_get_attribute(data,att[i]),_get_attribute(syn, att[i]))

			ax1.plot(log_tau, std, label=r"$\sigma$" + labels[i], color='#FF2C00')
			
			# Set limits
			ax1.set_xlim(log_tau[0], lim_max)

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{500}$")
			ax1.set_ylabel(r"$\sigma$" + labels_y[i])
			
			plt.tight_layout()

			plt.savefig(savepath + "analysis_" + inputs[i][1:] + add)
			
			if "-v" in sys.argv:
				print("Parameter:",inputs[i][1:])
				print("| log τ  | σ    |")
				print("-------------------")
				print("|  %1.3f | %.1f |" % (log_tau[np.abs(log_tau - ( 0)).argmin()],   std[np.abs(log_tau - ( 0)).argmin()]))
				print("| %1.3f | %.1f |" %  (log_tau[np.abs(log_tau - (-0.5)).argmin()], std[np.abs(log_tau - (-0.5)).argmin()]))
				print("| %1.3f | %.1f |" %  (log_tau[np.abs(log_tau - (-1)).argmin()],   std[np.abs(log_tau - (-1)).argmin()]))
				print("| %1.3f | %.1f |" %  (log_tau[np.abs(log_tau - (-1.5)).argmin()], std[np.abs(log_tau - (-1.5)).argmin()]))
				print("| %1.3f | %.1f |" %  (log_tau[np.abs(log_tau - (-2)).argmin()],   std[np.abs(log_tau - (-2)).argmin()]))
				temp = np.argmin(std)
				print(f"Minima at {log_tau[temp]} with {std[temp]}")
				print()
				print()


	if True:
		if "-vertical" in sys.argv:
			fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16), sharex = True,
														gridspec_kw = dict(hspace=0))
			fig.subplots_adjust(hspace=0, wspace=0)
		else:
			fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

		# Plot uncertainties
		# Standard deviation
		stdT = _std_dev(data.T[:,0,:], syn.T[:,0,:])
		ax1.plot(log_tau, stdT, label=r"$\sigma$" + labels_y[1], color = colors[0])

		stdB = _std_dev(data.B[:,0,:], syn.B[:,0,:])
		ax2.plot(log_tau, stdB, label=r"$\sigma$" + labels_y[4], color = colors[1])

		std = _std_dev(data.vlos[:,0,:], syn.vlos[:,0,:])
		ax3.plot(log_tau, std, label=r"$\sigma$" + labels_y[5], color = colors[2])

		std = _std_dev(data.gamma[:,0,:], syn.gamma[:,0,:])
		ax4.plot(log_tau, std, label=r"$\sigma$" + labels_y[6], color = colors[3])

		# Maximum difference
		if "-max" in sys.argv:
			ax1.plot(log_tau, np.min(data.T[:,0,:]-syn.T[:,0,:], axis = 0), "--",
					 label = r"$\Delta $" + labels_y[1] + r"$^{\mathrm{max}}$", color=colors[0])
			ax1.plot(log_tau, np.max(data.T[:,0,:]-syn.T[:,0,:], axis = 0), "--", color=colors[0])
			ax2.plot(log_tau, np.min(data.B[:,0,:]-syn.B[:,0,:], axis = 0), "--",
					 label = r"$\Delta $" + labels_y[4] + r"$^{\mathrm{max}}$", color=colors[1])
			ax2.plot(log_tau, np.max(data.B[:,0,:]-syn.B[:,0,:], axis = 0), "--", color=colors[1])
			ax3.plot(log_tau, np.min(data.vlos[:,0,:]-syn.vlos[:,0,:], axis = 0), "--",
					 label = r"$\Delta $" + labels_y[5] + r"$^{\mathrm{max}}$", color=colors[2])
			ax3.plot(log_tau, np.max(data.vlos[:,0,:]-syn.vlos[:,0,:], axis = 0), "--", color=colors[2])
			ax4.plot(log_tau, np.min(data.gamma[:,0,:]-syn.gamma[:,0,:], axis = 0), "--",
					 label = r"$\Delta $" + labels_y[6] + r"$^{\mathrm{max}}$", color=colors[3])
			ax4.plot(log_tau, np.max(data.gamma[:,0,:]-syn.gamma[:,0,:], axis = 0), "--", color=colors[3])

		# Set limits
		ax1.set_xlim(log_tau[0], lim_max)
		ax2.set_xlim(log_tau[0], lim_max)
		ax3.set_xlim(log_tau[0], lim_max)
		ax4.set_xlim(log_tau[0], lim_max)

		# Set labels
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r"$\log \tau_{c}$")
		ax4.set_xlabel(r"$\log \tau_{c}$")

		ax1.set_ylabel(r"$\sigma$" + labels[1])
		ax2.set_ylabel(r"$\sigma$" + labels[4])
		ax3.set_ylabel(r"$\sigma$" + labels[5])
		ax4.set_ylabel(r"$\sigma$" + labels[6])

		if "-vertical" not in sys.argv:
			ax1.set_title(titles[1])
			ax2.set_title(titles[4])
			ax3.set_title(titles[5])
			ax4.set_title(titles[6])
		
		if "-vertical" in sys.argv:
			xtitle = 0.41
		else:
			xtitle = 0.55
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle)

		if "-vertical" in sys.argv:	
			plt.tight_layout(pad=2.5,h_pad=0.0)
		else:
			plt.tight_layout(pad=1.5)

		plt.savefig(savepath + "analysis_result" + add)

	#######################################################################
	#					 Print out information				  #
	#######################################################################
		if "-v" in sys.argv:
			print("| log τ  | σ_T\t\t| σ_B \t  |")
			print("------------------------------------")
			print("|  %1.3f | %2.1f K\t| %.1f G |" % (log_tau[np.abs(log_tau - ( 0)).argmin()],   stdT[np.abs(log_tau - ( 0)).argmin()],   stdB[np.abs(log_tau - ( 0)).argmin()]))
			print("| %1.3f | %2.1f K\t| %.1f G |" % (log_tau[np.abs(log_tau - (-0.5)).argmin()], stdT[np.abs(log_tau - (-0.5)).argmin()], stdB[np.abs(log_tau - (-0.5)).argmin()]))
			print("| %1.3f | %2.1f K\t| %.1f G |" % (log_tau[np.abs(log_tau - (-1)).argmin()],   stdT[np.abs(log_tau - (-1)).argmin()],   stdB[np.abs(log_tau - (-1)).argmin()]))
			print("| %1.3f | %2.1f K\t| %.1f G |" % (log_tau[np.abs(log_tau - (-1.5)).argmin()], stdT[np.abs(log_tau - (-1.5)).argmin()], stdB[np.abs(log_tau - (-1.5)).argmin()]))
			print("| %1.3f | %2.1f K\t| %.1f G |" % (log_tau[np.abs(log_tau - (-2)).argmin()],   stdT[np.abs(log_tau - (-2)).argmin()],   stdB[np.abs(log_tau - (-2)).argmin()]))
			temp = np.argmin(stdB)
			print(f"Minima at {log_tau[temp]} with {stdB[temp]} G")


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
	conf = sir.read_config(sys.argv[1])
	if conf["mode"] == "MC":
		analysis(conf)
	else:
		print(f"[visualizer] Mode '{conf['mode']}' unknown or not defined")
	
		

