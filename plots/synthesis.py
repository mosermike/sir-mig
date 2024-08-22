"""
Plots a single inversion
"""

import numpy as np 
import sys
import os
from os.path import exists
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../src")) 
import sir
import definitions as d
import model_atm as m
import profile_stk as p

def _help():
	"""
	Help Page
	"""
	print("synthesis - Plots the result of a synthesis")
	print("Usage: python synthesis.py [OPTION]")
	print()
	sir.option("[1. Pos]","Config File")
	sir.option("[2. Pos]","Model number")
	
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-add_label","Additional label in plots")
	sir.option("-title","Title in Stokes plot")
	sir.option("-T","Plot temperature in K")
	sir.option("-Pe","Plot electron pressure in dyn/cm^2")
	sir.option("-vmicro","Plot microturbulence in cm/s")
	sir.option("-B","Plot magentic field strength in Gauss")
	sir.option("-vlos","Plot line of sight velocity in cm/s")
	sir.option("-gamma","Plot inclination by subtracting in deg")
	sir.option("-phi","Plot azimuth by adding in deg")
	sir.option("-z","Plot real height in km")
	sir.option("-Pg","Plot gas pressure in dyn/cm^2")
	sir.option("-rho","Plot density")
	sir.option("-syn","Synthesised model .mod file")
	sir.option("-vertical","Plot spectra vertically")
	sys.exit()

def synthesis(conf : dict, num : int):
	""""
	Plots the result of one inversion

	Parameters
	----------
	conf : dict
		Config. infos
	num : int
		Model number

	
	Returns
	-------
	None

	Other Parameters
	----------------
	Additional parameters given as an argument when the script is executed.
	-T
		Plot temperature in K
	-Pe
		Plot electron pressure in dyn/cm^2
	-vmicro
		Plot microturbulence in cm/s
	-B
		Plot magentic field strength in Gauss
	-vlos
		Plot line of sight velocity in km/s
	-gamma
		Plot inclination by subtracting in deg
	-phi
		Plot azimuth by adding in deg
	-z
		Plot real height in km
	-Pg
		Plot gas pressure in dyn/cm^2
	-rho
		Plot density
	-vertical
		Plot spectra vertically
	-save [str]
		Savepath, if not actual directory desired
	-add [str]
		Additional string at the end of file names
	-title [str]
		Title of the figures
	-add_label [str]
		Additional label in the legend. Can be used e.g. as a title
	-vertical
		Plot it vertically
	-num [int]
		which line number is used (only for mode `MC`)

	Raises
	------
	IndexError
		if num is out of range
	NotImplementedError
		if mode is not 'SY' or 'MC'
	"""
	# Import library
	sir.mpl_library()

	
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	# Check range
	if conf["mode"] != "SY" or conf["mode"] == "MC":
		raise NotImplementedError("[synthesis] Mode not defined. Select 'MC' or 'SY'!")
	

	path1 = conf["path"]
	if conf['mode'] == "MC":
		syn1 = m.read_model(os.path.join(path1,conf["syn_out"] + d.end_models))	# Synthesis Models 1

	elif conf['mode'] == "SY":
		syn1 = m.read_model(os.path.join(path1,conf["syn_in"]))	# Synthesis Models 1

	if conf['mode'] == "MC":
		obs1 = p.read_profile(os.path.join(path1,conf["syn_out"] + d.end_stokes))	# Synthesis Profiles 1
	elif conf['mode'] == "SY":
		obs1 = p.read_profile(os.path.join(path1,conf["syn_out"]))	# Synthesis Profiles 1
	
	if(num > syn1.nx):
		raise IndexError(f"[synthesis] num {num} out of range!")
	
	# Write in absolute wavelengths
	obs1.transform_wave_sir_to_abs(os.path.join(conf["path"],conf["lines"]))
	obs1.data_cut_wave = True

	# Observation from synthesis
	ll1, I1, Q1, U1, V1 = obs1.wave, obs1.stki[num,0],obs1.stkq[num,0],obs1.stku[num,0],obs1.stkv[num,0]
	
	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = sys.argv[sys.argv.index("-save")+1]
		if "/" in savepath:
			if not exists(savepath[:savepath.rfind("/")]):
				os.mkdir(savepath[:savepath.rfind("/")])
	
	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Title
	title = ''
	if '-title' in sys.argv:
		title = sys.argv[sys.argv.index("-title")+1]
	
	# Add label
	add_label = '_'
	if '-add_label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-add_label")+1]


	label_x = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if label_x != '0':
		ll1 -= float(label_x)
		ll2 -= float(label_x)

	########################
	#  Plot I, Q, U and V  #
	########################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,14), sharex=True,
			 gridspec_kw=dict(hspace=0), layout="compressed")
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), layout="compressed")


	# Add label
	ax2.scatter([], [], color="w", alpha=0, label=add_label)

	############################
	# Plot the Stokes profiles #
	############################
	ax1.plot(ll1, I1, "-", label=llabel)
	ax2.plot(ll1, Q1, "-", label=llabel)
	ax3.plot(ll1, U1, "-", label=llabel)
	ax4.plot(ll1, V1, "-", label=llabel)

	# Set xlimits
	ax1.set_xlim(ll1[0], ll1[-1])
	ax2.set_xlim(ll1[0], ll1[-1])
	ax3.set_xlim(ll1[0], ll1[-1])
	ax4.set_xlim(ll1[0], ll1[-1])

	ax1.set_ylabel(r'$I / I_c$')
	ax2.set_ylabel(r'$Q / I_c$')
	ax3.set_ylabel(r'$U / I_c$')
	ax4.set_ylabel(r'$V / I_c$')
	
	if label_x != "0":
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' [\AA]', loc='center')
	else:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\lambda$' + r' [\AA]', loc='center')
		ax4.set_xlabel(r'$\lambda$' + r' [\AA]', loc='center')
	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle1 = 0.41
		else:
			xtitle1 = 0.5
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle1)
		elif title is not None:
			fig.suptitle(title, y=0.98, x=xtitle1)

	#########################
	# Set Legend and Limits #
	#########################
	if "-vertical" in sys.argv:
		if (conf["mode"] != "SY"):
			ax1.set_ylim(0.9*np.min(np.abs(I1)) , 1.1*np.max(np.abs(I1)))
			ax2.set_ylim(-1.1*np.max(np.abs(Q1)), 1.1*np.max(np.abs(Q1)))
			ax3.set_ylim(-1.1*np.max(np.abs(U1)), 1.1*np.max(np.abs(U1)))
			ax4.set_ylim(-1.1*np.max(np.abs(V1)), 1.1*np.max(np.abs(V1)))
		#ax1.legend(bbox_to_anchor=(1.01,0.95))
		ax1.legend()
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))

	
	plt.savefig(savepath + "synthesis_stokes_" + str(num) + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-gamma", "-phi", "-z", "-Pg","-rho"]
	index  = [1,2,3,4,5,6,7,8,9,10]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$"]
	titles   = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength B", r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
			 r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$"]

	i = 0
	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Get colors used in the actual cycle
	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8,7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)

			# Plot
			# Synthesised model = Real model
			ax1.plot(syn1.tau, syn1.get_attribute(inputs[i][1:])[num,0], label=f"{llabel}",color=colors[0])

			# Set xlimits
			ax1.set_xlim(syn1.tau[0], syn1.tau[-1])

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{c}$")
			ax1.set_ylabel(labels[index[i]])

			if index[i] == 2 or index[i] == 9:
				ax1.semilogy()

			ax1.set_title(titles[i])
			# set the spacing between subplots
			fig.savefig(savepath + "synthesis_model_" + str(num) + "_" + str(inputs[i][1:]) + add)
		
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0), layout="compressed")
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), layout="compressed")

	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Get colors used in the actual cycle
	llabel = "Syn. Model"
	ax1.plot(syn1.tau, syn1.T[num,0], label=f"{llabel}", color=colors[0])
	ax2.plot(syn1.tau, syn1.B[num,0], label=f"{llabel}", color=colors[0])
	ax3.plot(syn1.tau, syn1.vlos[num,0], label=f"{llabel}", color=colors[0])
	ax4.plot(syn1.tau, syn1.gamma[num,0], label=f"{llabel}", color=colors[0])
	
	#####################
	#	Set limits	#
	#####################
	ax1.set_xlim(syn1.tau[0], syn1.tau[-1])
	ax2.set_xlim(syn1.tau[0], syn1.tau[-1])
	ax3.set_xlim(syn1.tau[0], syn1.tau[-1])
	ax4.set_xlim(syn1.tau[0], syn1.tau[-1])

	Min, Max = ax2.get_ylim()
	ax2.set_ylim(Min,Max*1.15)
	Min, Max = ax3.get_ylim()
	ax3.set_ylim(Min,Max*1.15)
	Min, Max = ax4.get_ylim()
	ax4.set_ylim(Min,Max*1.15)

	#####################
	#	Set labels	#
	#####################
	if "-vertical" not in sys.argv:
		ax1.set_xlabel(r"$\log \tau_{c}$")
		ax2.set_xlabel(r"$\log \tau_{c}$")
		ax3.set_xlabel(r"$\log \tau_{c}$")
	ax4.set_xlabel(r"$\log \tau_{c}$")

	ax1.set_ylabel(labels[1])
	ax2.set_ylabel(labels[4])
	ax3.set_ylabel(labels[5])
	ax4.set_ylabel(labels[6])

	###############################
	#	Set title and legend	#
	###############################
	if "-vertical" not in sys.argv:
		ax1.set_title(titles[0])
		ax2.set_title(titles[3])
		ax3.set_title(titles[4])
		ax4.set_title(titles[5])

	# Set title position depending on the chosen plot and consider the flags hinode and gris
	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle2 = 0.51
		else:
			xtitle2 = 0.55
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle2)

	fig.savefig(savepath + "synthesis_result_" + str(num) + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
	conf = sir.read_config(sys.argv[1])
	synthesis(conf, int(sys.argv[2]),int(sys.argv[3]))




