"""

Analyses the quality of multiple Monte Carlo Simulations

"""

import numpy as np 
import sys
import os
sys.path.append(sys.path[0] + "/../")
sys.path.append(sys.path[0] + "/../tools")
import sir
import definitions as d
from model_atm import *  # Model class
import matplotlib.pyplot as plt
from os.path import exists
import os

#############
# Help page #
#############
def _help():
	"""
	Help Page
	"""
	print("analysis_multiple - Analyses the quality of multiple Monte Carlo Simulations")
	print("Usage: python analysis_multiple.py [OPTION]")
	print()
	sir.option("[1. Pos]","Config Files as a list config1,config2,config3,...")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-T","Compare temperature in K")
	sir.option("-B","Compare magentic field strength in Gauss")
	sir.option("-vlos","Compare line of sight velocity in km/s")
	sir.option("-inc","Compare inclination by subtracting in deg")
	sir.option("-azi","Compare azimuth by adding in deg")
	sir.option("-title","Title of the 4 figures plot")
	sir.option("-xtitle","x position of title in Stokes plot (optional)")
	sir.option("-limitT","Set y limits in T as 'ymin,ymax'")
	sir.option("-vertical","Plot the last plot vertical")
	sir.option("-v","print out tables with values at different log taus.")
	print()
	print("Note: B, vlos, inc and T is always compared but not plotted alone if the flags are not used.")
	sys.exit()


def _get_attribute(model, att):
	'''
	Returns a model parameter depending on a string

	Parameter
	---------
	model : class Model
	att : str

	'''
	if att == "tau":
		return model.tau
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
	
	
def analysis_multiple(confs, labels):
	"""
	Analysis multiple simulations and puts it in one plot

	Parameters
	----------
	confs : list
		List of dictionaries with the configs
	labels : list
		List of strings with the labels for each plot

	Returns
	-------
	None

	Other Parameters
	----------------
	There are optional options which change the plots. Additional parameters given as an argument when the script is executed.

	-save [str], optional
		Additional save path. Default './'
	-add [str]
		Additional text in filenames
	-T
		Compare temperature in K
	-B
		Compare magentic field strength in Gauss
	-vlos
		Compare line of sight velocity in km/s
	-inc
		Compare inclination by subtracting in deg
	-azi
		Compare azimuth by adding in deg
	-title [str]
		Title of the 4 figures plot
	-xtitle [float]
		x position of title in Stokes plot
	-limitT
		Set y limits in T as 'ymin,ymax'
	-vertical
		Plot the last plot vertical
	-v
		print out tables with values at different log taus.
	
	Notes
	-----
	Note: B, vlos, inc and T is always compared but not plotted alone if the flags are not used.

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


	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]  # Get colors used in the actual cycle
	
	limitT = (None,None)
	if "-limitT" in sys.argv:
		limitT = [float(i) for i in sys.argv[sys.argv.index("-limitT")+1].split(",")]

	linestyle_str = [
		'solid',      # Same as (0, ()) or '-'
		'dotted',    # Same as (0, (1, 1)) or ':'
		'dashed',    # Same as '--'
		'dashdot',
		'solid',      # Same as (0, ()) or '-'
		'dotted',    # Same as (0, (1, 1)) or ':'
		'dashed',    # Same as '--'
		'dashdot',
		(0, (3,10,1,10))
	]

	#######################################################################
	#	    READ INPUT, WRITE DEFINITIONS AND LOAD DATA			  	#
	#######################################################################

	# Additional labels
	save = ''
	if "-save" in sys.argv:
		save = sys.argv[sys.argv.index("-save")+1]
		# Create savepath
		if not exists(save[0:save.rfind('/')]):
			os.mkdir(save[0:save.rfind('/')])

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

	#######################################################################
	#		    DETERMINE WHAT PARAMETER ARE CHOSEN				#
	#######################################################################
	inputs = ["__", "-T", "-Pe", "-vmicro", "-B", "-vlos", "-inc", "-azi", "-z", "-Pg", "-rho"]
	att = ["tau", "T", "Pe", "vmicro", "B", "vlos", "gamma", "phi", "z", "Pg", "rho"]  # For getting the value from the class
	labels_y = ["__", r"$_T$ [K]", "", "", r"$_B$ [G]",
					r"$_{\mathrm{v}_{\mathrm{los}}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$",
					r"$_\gamma$ [deg]", r"$_\phi$ [deg]", r"$_z$ [km]",
					r"$_{P_g}$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$",
					r"$_{rho}$ $\left[\frac{\mathrm{g}}{\mathrm{cm}^3}\right]"
				]
	Labels = ["__", r"$_T$ [K]", "", "", r"$_B$ [G]",
				r"$_{\mathrm{v}_{\mathrm{los}}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$",
				r"$_\gamma$ [deg]", r"$_\phi$ [deg]", r"$_z$ [km]",
				r"$_{P_g}$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$",
				r"$_{\rho}$ $\left[\frac{\mathrm{g}}{\mathrm{cm}^3}\right]$"
			]
	titles = ["", r"Temperature T", "", "", r"Magnetic Field Strength B",
				r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$", r"Inclination $\gamma$",
				r"Azimuth $\phi$", r"Height z", r"Gas Pressure $P_g$"
			]

	############################
	# Load physical parameters #
	############################
	# Create empty arrays for the differences
	log_taus = []
	syns = []
	fits = []
	nums = []

	# Go through each model and load models
	for n in range(len(confs)):
		path = confs[n]["path"]
		num = confs[n]['num']
		fit = read_model(os.path.join(path, confs[n]["inv_out"]) + d.end_models)
		syn = read_model(os.path.join(path, confs[n]["syn_out"]+ d.end_models))

		# Correct phi range
		syn.correct_phi()
		fit.correct_phi()

		log_tau = fit.tau

		log_taus.append(log_tau)
		syns.append(syn)
		fits.append(fit)
		nums.append(num)

	############################
	# Plot physical parameters #
	############################

	# Determine maximum of plot in x and cut data if needed
	lim_max = -3
	for i in range(len(fits)):
		fits[i].set_limit(lim_max)
		syns[i].set_limit(lim_max)
	lim_max = fits[0].tau[-1]
	#####
	# Plot the physical parameters separately as mentioned in the flags
	#####

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			stds = []
			if len(labels) > 5:
				fig, ax1 = plt.subplots(figsize=(9, 7), layout="compressed")
			else:
				fig, ax1 = plt.subplots(figsize=(8, 7), layout="compressed")

			plt.title(titles[i])

			for n in range(len(log_taus)):
				# Standard deviation
				std = np.sqrt(np.sum((_get_attribute(fits[n],att[i]) - _get_attribute(syns[n], att[i]))**2, axis=0)/(nums[n]-1))	 
				ax1.plot(fits[n].tau, std, label=labels[n], linestyle=linestyle_str[n % len(linestyle_str)])
				stds.append(std)
			# Set limits
			ax1.set_xlim(fits[0].tau[0], lim_max)

			if inputs[i] == "-T":
				ax1.set_ylim(limitT)

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{c}$")
			ax1.set_ylabel(r"$\sigma$" + labels_y[i])
			
			if len(labels) > 5:
				ax1.legend(loc='center right', bbox_to_anchor=(1.25+np.max([len(i) for i in labels])/85, 0.5), frameon=False)
			else:
				ax1.legend()
			#plt.tight_layout()

			fig.savefig(save + "analysis_multiple_" + inputs[i][1:] + add)

			if "-v" in sys.argv:
				print("Parameter:", inputs[i][1:])
				print("| log τ  |", end="")
				for j in range(len(fits)):
					print(f" σ_{labels[j]}    |", end="")
				print()
				print("-----------------------")
				for i in [1, 0.5, 0.0, -0.5, -1, -1.5, -2]:
					print("|  %1.2f |" % (fits[0].tau[np.abs(fits[0].tau - (i)).argmin()]), end="")
					for j in range(len(fits)):
						print("  %1.2f |" % (stds[j][np.abs(fits[j].tau[0] - (i)).argmin()]), end="")
					print()
				print()
				for j in range(len(fits)):
					temp1 = np.argmin(stds[j])
					print(f"Minima 1 at {fits[j].tau[0][temp1]} with {'%1.2f' % (stds[j][temp1])}")
				print()
				print()		


	###############################
	#	Plot uncertainties		#
	###############################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots( 4, 1, figsize=(12,16), sharex=True,
												gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if len(labels) > 5:
			fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
		else:
			fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))

	######################
	# Standard deviation #
	######################
	stdTs = []
	stdBs = []
	for n in range(len(log_taus)):	
		std = np.sqrt(np.sum((fits[n].T[:,0,:] - syns[n].T[:,0,:])**2, axis=0)/(nums[n]-1))
		stdTs.append(std)
		ax1.plot(fits[n].tau, stdTs[n], label=labels[n], linestyle=linestyle_str[n % len(linestyle_str)])

		std = np.sqrt(np.sum((fits[n].B[:,0,:] - syns[n].B[:,0,:])**2, axis=0)/(nums[n]-1))
		stdBs.append(std)
		ax2.plot(fits[n].tau, stdBs[n], label=labels[n], linestyle=linestyle_str[n % len(linestyle_str)])

		std = np.sqrt(np.sum((fits[n].vlos[:,0,:] - syns[n].vlos[:,0,:])**2, axis=0)/(nums[n]-1))
		ax3.plot(fits[n].tau, std, label=labels[n], linestyle=linestyle_str[n % len(linestyle_str)])

		std = np.sqrt(np.sum((fits[n].gamma[:,0,:] - syns[n].gamma[:,0,:])**2, axis=0)/(nums[n]-1))
		ax4.plot(fits[n].tau, std, label=labels[n], linestyle=linestyle_str[n % len(linestyle_str)])

	##############
	# Set limits #
	##############
	ax1.set_xlim(fits[0].tau[0], lim_max)
	ax2.set_xlim(fits[0].tau[0], lim_max)
	ax3.set_xlim(fits[0].tau[0], lim_max)
	ax4.set_xlim(fits[0].tau[0], lim_max)

	if "-limitT" in sys.argv:
		ax1.set_ylim(limitT)

	##############
	# Set labels #
	##############
	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r"$\log \tau_{c}$")
	ax4.set_xlabel(r"$\log \tau_{c}$")
	ax1.set_ylabel(r"$\sigma$" + Labels[1])
	ax2.set_ylabel(r"$\sigma$" + Labels[4])
	ax3.set_ylabel(r"$\sigma$" + Labels[5])
	ax4.set_ylabel(r"$\sigma$" + Labels[6])

	#########
	# Title #
	#########
	# Set axis titles if not in vertical
	if "-vertical" not in sys.argv:
		ax1.set_title(titles[1], fontsize=20)
		ax2.set_title(titles[4], fontsize=20)
		ax3.set_title(titles[5], fontsize=20)
		ax4.set_title(titles[6], fontsize=20)

	###############
	# Plot legend #
	###############
	if len(labels) > 5:
		ax2.legend(loc='center right', bbox_to_anchor=(1.4 + len(labels[0])/100, 0.5), frameon=False)
	else:
		ax1.legend()

	#########################
	# Title and positioning #
	#########################
	# Place the title depending on the wanted design of the plot
	if "-vertical" in sys.argv:
		xtitle = 0.41
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle)
	else:
		if len(title) < 50:
			xtitle = 0.55
		else:
			xtitle = 0.45
		xtitle = 0.5
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle)	    

	if "-vertical" in sys.argv:	
		plt.tight_layout(pad=2.5,h_pad=0.0)
	else:
		plt.tight_layout(pad=1.5)

	plt.savefig(save + "analysis_multiple" + add)

	#########################
	# Print out information #
	#########################
	if "-v" in sys.argv:
		print("| log τ  |", end="")

		for n in range(len(fits)):
			print(f" σ_T_{labels[n]} in K\t|", end="")
		for n in range(len(fits)):
			print(f" σ_B_{labels[n]} in G\t|", end="")
		print()
		print("--------------------------------------------------------------------------------------------")
		for i in [1, 0.5, 0, -0.5, -1, -1.5, -2]:
			print("|  %1.2f \t\t |" % (fits[0].tau[np.abs(fits[0].tau - i).argmin()]), end="")
			for j in range(len(labels)):
				for n in range(len(fits)):
					if j == 0:
						print("  %1.1f \t\t |" % (stdTs[n][np.abs(fits[n].tau - i).argmin()]), end="")
					else:
						print("  %1.1f \t\t |" % (stdBs[n][np.abs(fits[n].tau - i).argmin()]), end="")
			print()

		print()

		for j in range(len(labels)):
			temp1 = np.argmin(stdBs[j])
			print(f"Minima {j+1} at {'%.2f' %fits[j].tau[temp1]} with {'%.3f' % stdBs[j][temp1]} G")

	return
##################################################################


if __name__ == "__main__":
	# Print help page
	if "-h" in sys.argv:
		_help()

	# Config files
	if "," in sys.argv[1]:
		confs = sys.argv[1].replace(", ", ",").split(',')
		confs = [sir.read_config(i, check=False) for i in confs]
		temp = sys.argv[1].split(",")
		for i in range(len(temp)):
			# Correct for ./ paths
			if confs[i]['path'] == "./":
				confs[i]['path'] = temp[i][:temp[i].rfind('/')+1]
	else:
		confs  = [sir.read_config(sys.argv[1])]

	
	if confs[0]["mode"] == "MC":
		#Labels
		labels = []
		for i in range(len(confs)):
			ll = input(f"Label {i+1}: ")
			labels.append(ll)
		analysis_multiple(confs,labels)
	else:
		print(f"[visualizer] Mode '{confs[0]['mode']}' unknown or not defined")

