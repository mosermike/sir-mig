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
	print("inversion - Plots the result of one inversion")
	print("Usage: python inversion.py [OPTION]")
	print()
	sir.option("[1. Pos]","Config File")
	sir.option("[2. Pos]","x position")
	sir.option("[3. Pos]","y position (choose '0' for mode 'MC')")
	
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
	sir.option("-num:","Number of the line considered (Default: take first one) (for Mode 'MC')")
	sys.exit()

def inversion(conf : dict, x : int, y : int):
	""""
	Plots the result of one inversion

	Parameters
	----------
	conf : dict
		Config. infos
	x : int
		Absolute x position inversion
	y : int
		Absolute y position inversion

	
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
	"""
	# Import library
	sir.mpl_library()

	
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################

	path1 = conf["path"]
	if conf['mode'] == "MC":
		syn1 = m.read_model(os.path.join(path1,conf["syn_out"] + d.end_models))	# Synthesis Models 1

	if conf1['mode'] == "2C":
		syn1 = m.read_model(os.path.join(path1,conf["inv_out"] + d.end_models1))	# Synthesis Models 2	
		err11 = m.read_model(os.path.join(path1,conf["inv_out"] + d.end_errors1))	# Inversion Models 2
		phy1 = m.read_model(os.path.join(path1,conf["inv_out"] + d.end_models2))	# Inversion Models 2
		err1 = m.read_model(os.path.join(path1,conf["inv_out"] + d.end_errors2))	# Inversion Models 2
	else:
		phy1 = m.read_model(os.path.join(path1,conf["inv_out"] + d.end_models))	# Inversion Models 1
		err1 = m.read_model(os.path.join(path1,conf["inv_out"] + d.end_errors))	# Inversion Models 1
		
	if conf1['mode'] == "MC":
		obs1 = p.read_profile(os.path.join(path1,conf["syn_out"] + d.end_stokes))	# Synthesis Profiles 1
	else:
		obs1 = p.read_profile(os.path.join(path1,conf["cube"]))
	fit1 = p.read_profile(os.path.join(path1,conf["inv_out"] + d.end_stokes))	# Inversion Profiles 1



	# Cut wave
	if conf['mode'] == "1C" or conf["mode"] == "2C":
		obs1.cut_to_wave(conf["range_wave"])
		obs1.cut_to_map(conf["map"])

	else:
		num = phy1.indx[0]
		if "-num" in sys.argv:
			num = int(sys.argv[sys.argv.index("-num")+1])
		# Cut the wave to the line number
		ll1 = np.where(phy1.indx==num)[0][0]
		ll2 = np.where(phy1.indx==num)[0][-1]
		obs1.cut_to_wave([obs1.wave[ll1],obs1.wave[ll2]])
		fit1.cut_to_wave([fit1.wave[ll1],fit1.wave[ll2]])

	# Change to relative x and y for mode 1C and 2C
	if conf["mode"] == "1C" or conf["mode"] == "2C":
		x = x - conf["map"][0]
		y = y - conf["map"][2]

	# Observation from synthesis
	ll1, I1, Q1, U1, V1 = obs1.wave, obs1.stki[x,y],obs1.stkq[x,y],obs1.stku[x,y],obs1.stkv[x,y]
	
	# Load fit data
	ll2, I_fit1, Q_fit1, U_fit1, V_fit1 = fit1.wave, fit1.stki[x,y],fit1.stkq[x,y],fit1.stku[x,y],fit1.stkv[x,y]

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath):
			os.mkdir(savepath)
	
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



	#############################################################
	#					PLOTTING STUFF					#
	#############################################################
	# Change to abs. wavelength to the line core of the first number
	if conf1['mode'] == "MC":
		ll0 = input("Put wavelength of the line core: ")
		ll1 += ll0
		ll2 += ll0

	label_x = "0"

	temp = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if temp != '0':
		ll1 -= float(temp)
		ll2 -= float(temp)
		label_x = temp
	elif conf["mode"] == "MC":
		# Keep the same wavelengths as in the Grid file
		label_x = str(ll0)
		ll1 -= ll0
		ll2 -= ll0


	########################
	#  Plot I, Q, U and V  #
	########################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,14), sharex=True,
			 gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12))


	# Add label
	ax2.scatter([], [], color="w", alpha=0, label=add_label)

	############################
	# Plot the Stokes profiles #
	############################
	if conf1["mode"] == "MC":
		llabel = "Syn. Profile"
	else:
		llabel = "Obs. Profile"
	ax1.plot(ll1, I1, "x", label=llabel)
	ax2.plot(ll1, Q1, "x", label=llabel)
	ax3.plot(ll1, U1, "x", label=llabel)
	ax4.plot(ll1, V1, "x", label=llabel)

	ax1.plot(ll1, I_fit1, "-", label = "Best Fit")
	ax2.plot(ll1, Q_fit1, "-", label = "Best Fit")
	ax3.plot(ll1, U_fit1, "-", label = "Best Fit")
	ax4.plot(ll1, V_fit1, "-", label = "Best Fit")

	# Set xlimits
	ax1.set_xlim(ll1[0], ll1[-1])
	ax2.set_xlim(ll1[0], ll1[-1])
	ax3.set_xlim(ll1[0], ll1[-1])
	ax4.set_xlim(ll1[0], ll1[-1])

	#######################################################################
	# Set labels												#
	# The labels depend on the arguments and the chosen line			#
	# The code is so long as I tried to make it as automatic as possible	#
	#######################################################################
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
		ax1.set_ylim(0.9*np.min(np.abs(np.append(I1,I_fit1))), 1.1*np.max(np.abs(np.append(I1,I_fit1))))
		ax2.set_ylim(-1.1*np.max(np.abs(np.append(Q1,Q_fit1))), 1.1*np.max(np.abs(np.append(Q1,Q_fit1))))
		ax3.set_ylim(-1.1*np.max(np.abs(np.append(U1,U_fit1))), 1.1*np.max(np.abs(np.append(U1,U_fit1))))
		ax4.set_ylim(-1.1*np.max(np.abs(np.append(V1,V_fit1))), 1.1*np.max(np.abs(np.append(V1,V_fit1))))
		ax1.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout(pad=2.5)

	plt.savefig(savepath + "inversion2_stokes_x" + str(x) + "_y" + str(y) + add)

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

			# '0C5DA5', 'FF2C00', '00B945', 'FF9500', '845B97', '474747', 'C20078','054907'
			# Plot
			# Synthesised model = Real model
			if conf1["mode"] == "MC" or conf1['mode'] == "2C":
				if conf1['mode'] == "2C":
					llabel = "Best Fit M. 1"
				else:
					llabel = "Syn Model"
				ax1.plot(syn1.tau, syn1.get_attribute(inputs[i][1:])[x,y], label=f"{llabel}",color=colors[0])
			if conf1['mode'] == "2C":
				llabel = "Best Fit M. 2"
			else:
				llabel = "Best Fit"
			
			ax1.plot(fit1.tau, fit1.get_attribute(inputs[i][1:])[x,y], label=f"{llabel}",color = colors[1])

			if conf1['mode'] == "2C":
				# Error of fit
				ax1.fill_between(syn1.tau, syn1.get_attribute(inputs[i][1:])[x,y] - err11.get_attribute(inputs[i][1:])[x,y],
							syn1.get_attribute(inputs[i][1:])[x,y] + err11.get_attribute(inputs[i][1:])[x,y], alpha = 0.5,
							color=colors[0], lw=0)
			
			# Error of fit
			ax1.fill_between(fit1.tau, fit1.get_attribute(inputs[i][1:])[x,y] - err1.get_attribute(inputs[i][1:])[x,y],
						 fit1.get_attribute(inputs[i][1:])[x,y] + err1.get_attribute(inputs[i][1:])[x,y], alpha = 0.5,
						 color=colors[1], lw=0)



			# Set xlimits
			ax1.set_xlim(fit1.tau[0], fit1.tau[-1])

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{c}$")
			ax1.set_ylabel(labels[index[i]])

			if index[i] == 2 or index[i] == 9:
				ax1.semilogy()

			# Legend
			ax1.legend()
		
			ax1.set_title(titles[i])
			# set the spacing between subplots
			plt.tight_layout(pad=2)
			plt.savefig(savepath + "inversion_x" + str(x) + "_y"  + str(y)+ str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = phy1.tau[-1]
	if conf1["mode"] == "MC" or conf1['mode'] == "2C":
		syn1.set_limit(lim_max)
		syn1.set_limit(lim_max)
	phy1.set_limit(lim_max)
	err1.set_limit(lim_max)


	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	if conf1["mode"] == "MC":
		llabel = "Syn. Model"
		ax1.plot(syn1.tau, syn1.T[x,y], label=f"{llabel}", color='#0C5DA5')
		ax2.plot(syn1.tau, syn1.B[x,y], label=f"{llabel}", color='#0C5DA5')
		ax3.plot(syn1.tau, syn1.gamma[x,y], label=f"{llabel}", color='#0C5DA5')
		ax4.plot(syn1.tau, syn1.vlos[x,y], label=f"{llabel}", color='#0C5DA5')

	llabel = "Best Fit"
	ax1.plot(phy1.tau, phy1.T[x,y], label=f"{llabel}", color='#FF2C00')
	ax2.plot(phy1.tau, phy1.B[x,y], label=f"{llabel}", color='#FF2C00')
	ax3.plot(phy1.tau, phy1.vlos[x,y], label=f"{llabel}", color='#FF2C00')
	ax4.plot(phy1.tau, phy1.gamma[x,y], label=f"{llabel}", color='#FF2C00')


	if conf1['mode'] == "2C":
		ax1.fill_between(syn1.tau, syn1.T[x,y] - err11.T[x,y],
				 syn1.T[x,y] + err11.T[x,y], alpha = 0.5,
				 color='#0C5DA5', lw=0)
		ax2.fill_between(syn1.tau, syn1.B[x,y] - err11.B[x,y],
					syn1.B[x,y] + err11.B[x,y], alpha = 0.5,
					color='#0C5DA5', lw=0)
		ax3.fill_between(syn1.tau, syn1.vlos[x,y] - err11.vlos[x,y],
					syn1.vlos[x,y] + err11.vlos[x,y], alpha = 0.5,
					color='#0C5DA5', lw=0)
		ax4.fill_between(syn1.tau, syn1.gamma[x,y] - err11.gamma[x,y],
					syn1.gamma[x,y] + err11.gamma[x,y], alpha = 0.5,
					color='#0C5DA5', lw=0)
		


	ax1.fill_between(phy1.tau, phy1.T[x,y] - err1.T[x,y],
				 phy1.T[x,y] + err1.T[x,y], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax2.fill_between(phy1.tau, phy1.B[x,y] - err1.B[x,y],
				 phy1.B[x,y] + err1.B[x,y], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax3.fill_between(phy1.tau, phy1.vlos[x,y] - err1.vlos[x,y],
				 phy1.vlos[x,y] + err1.vlos[x,y], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax4.fill_between(phy1.tau, phy1.gamma[x,y] - err1.gamma[x,y],
				 phy1.gamma[x,y] + err1.gamma[x,y], alpha = 0.5,
				 color='#FF2C00', lw=0)

	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(phy1.tau[0], lim_max)
	ax2.set_xlim(phy1.tau[0], lim_max)
	ax3.set_xlim(phy1.tau[0], lim_max)
	ax4.set_xlim(phy1.tau[0], lim_max)

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

	ax1.legend(loc='upper right')
	ax2.legend(loc='upper right')
	ax3.legend(loc='upper right')
	ax4.legend(loc='upper right')

	# Set title position depending on the chosen plot and consider the flags hinode and gris
	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle2 = 0.51
		else:
			xtitle2 = 0.55
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle2)

	if "-vertical" in sys.argv:	
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		# set the spacing between subplots
		plt.tight_layout(pad=2)

	  
	plt.savefig(savepath + "inversion_result_x" + str(x) + "_y"  + str(y) + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
	conf1 = sir.read_config(sys.argv[1])
	inversion(conf1, int(sys.argv[2]),int(sys.argv[3]))




