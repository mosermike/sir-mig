"""
Plots two single inversions
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

def help():
	"""
	Help Page
	"""
	print("plot_inversion_2 - Plots the result of two inversions")
	print("Usage: python plot_inversion_2 [OPTION]")
	print()
	sir.option("[1. Pos]","Config Model 1")
	sir.option("[2. Pos]","x position Profile 1")
	sir.option("[3. Pos]","y position Profile 1")
	sir.option("[4. Pos]","Config Model 2")
	sir.option("[5. Pos]","x position Profile 2")
	sir.option("[6. Pos]","y position Profile 2")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-add_label","Additional label in plots")
	sir.option("-label1","Label for 1st inversion")
	sir.option("-label2","Label for 2nd inversion")
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
	sir.option("-num","Number of the line considered (Default: take first one) (for Mode 'MC')")
	sir.option("-one","Plot observations only once")
	sir.option("-err","Plot errorbars")
	sys.exit()

def inversion_2(conf1 : dict, x1 : int, y1 : int, conf2 : dict, x2 : int, y2 : int):
	""""
	Plots the result of one inversion

	Parameters
	----------
	conf1 : dict
		Config infos for 1st inversion
	x1 : int
		x position of 1st inversion
	y1 : int
		y position of 1st inversion
	conf2 : dict
		Config infos for 2nd inversion
	x2 : int
		x position of 2nd inversion
	y2 : int
		y position of 2nd inversion
	
	Returns
	-------
	None

	Other Parameters
	----------------
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
	-label1 [str]
		Label for the 1st inversion. Default: "1st Inversion"
	-label2 [str]
		Label for the 2nd inversion. Default: "2nd Inversion"
	-vertical
		Plot it vertically
	-num [int]
		which line number is used (only for mode `MC`)
	-one
		Observations are the same and therefore only plotted once
	-err
		Plot errorbars
	"""
	# Import library
	sir.mpl_library()

	
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################

	path1 = conf1["path"]
	if conf1['mode'] == "MC":
		syn1 = m.read_model(os.path.join(path1,conf1["syn_out"] + d.end_models))	# Synthesis Models 1

	if conf1['mode'] == "2C":
		syn1 = m.read_model(os.path.join(path1,conf1["inv_out"] + d.end_models1))	# Synthesis Models 2	
		err11 = m.read_model(os.path.join(path1,conf1["inv_out"] + d.end_errors1))	# Inversion Models 2
		phy1 = m.read_model(os.path.join(path1,conf1["inv_out"] + d.end_models2))	# Inversion Models 2
		err1 = m.read_model(os.path.join(path1,conf1["inv_out"] + d.end_errors2))	# Inversion Models 2
	else:
		phy1 = m.read_model(os.path.join(path1,conf1["inv_out"] + d.end_models))	# Inversion Models 1
		err1 = m.read_model(os.path.join(path1,conf1["inv_out"] + d.end_errors))	# Inversion Models 1
		
	if conf1['mode'] == "MC":
		obs1 = p.read_profile(os.path.join(path1,conf1["syn_out"] + d.end_stokes))	# Synthesis Profiles 1
	else:
		obs1 = p.read_profile(os.path.join(path1,conf1["cube"]))
	fit1 = p.read_profile(os.path.join(path1,conf1["inv_out"] + d.end_stokes))	# Inversion Profiles 1


	path2 = conf2["path"]
	if conf1['mode'] == "MC":
		syn2 = m.read_model(os.path.join(path2,conf2["syn_out"] + d.end_models))	# Synthesis Models 2	
	if conf1['mode'] == "2C":
		syn2 = m.read_model(os.path.join(path2,conf2["inv_out"] + d.end_models1))	# Synthesis Models 2	
		err22 = m.read_model(os.path.join(path2,conf2["inv_out"] + d.end_errors1))	# Inversion Models 2
		phy2 = m.read_model(os.path.join(path2,conf2["inv_out"] + d.end_models2))	# Inversion Models 2
		err2 = m.read_model(os.path.join(path2,conf2["inv_out"] + d.end_errors2))	# Inversion Models 2
	else:
		phy2 = m.read_model(os.path.join(path2,conf2["inv_out"] + d.end_models))	# Inversion Models 2
		err2 = m.read_model(os.path.join(path2,conf2["inv_out"] + d.end_errors))	# Inversion Models 2
	if conf1['mode'] == "MC":
		obs2 = p.read_profile(os.path.join(path2,conf2["syn_out"] + d.end_stokes))	# Synthesis Profiles 2
	else:
		obs2 = p.read_profile(os.path.join(path2,conf2["cube"]))
	fit2 = p.read_profile(os.path.join(path2,conf2["inv_out"] + d.end_stokes))	# Inversion Profiles 2

	# Cut wave
	if conf1['mode'] == "1C" or conf1["mode"] == "2C":
		obs1.cut_to_wave(conf1["range_wave"])
		obs1.cut_to_map(conf1["Map"])
		obs2.cut_to_wave(conf2["range_wave"])
		obs2.cut_to_map(conf2["Map"])
	else:
		num = fit1.indx[0]
		if "-num" in sys.argv:
			num = int(sys.argv[sys.argv.index("-num")+1])
		ind1 = 0
		ind2 = 0

		for i in range(len(conf1["atoms"])):
			if str(num) in conf1["atoms"][i]:
				ind1 = i
		for i in range(len(conf2["atoms"])):
			if str(num) in conf2["atoms"][i]:
				ind2 = i

		# define map as it is used later for saving
		conf1["map"] = [0,0,0,0]
		conf2["map"] = [0,0,0,0]
		
		obs1.cut_to_wave(np.array([conf1["range_wave"][ind1]]))
		obs2.cut_to_wave(np.array([conf2["range_wave"][ind2]]))
		fit1.cut_to_wave(np.array([conf1["range_wave"][ind1]]))
		fit2.cut_to_wave(np.array([conf2["range_wave"][ind2]]))

	# Observation from synthesis
	ll1, I1, Q1, U1, V1 = obs1.wave, obs1.stki[x1,y1],obs1.stkq[x1,y1],obs1.stku[x1,y1],obs1.stkv[x1,y1]
	ll2, I2, Q2, U2, V2 = obs2.wave, obs2.stki[x2,y2],obs2.stkq[x2,y2],obs2.stku[x2,y2],obs2.stkv[x2,y2]
	
	# Load fit data
	ll1_fit, I_fit1, Q_fit1, U_fit1, V_fit1 = fit1.wave, fit1.stki[x1,y1],fit1.stkq[x1,y1],fit1.stku[x1,y1],fit1.stkv[x1,y2]
	ll2_fit, I_fit2, Q_fit2, U_fit2, V_fit2 = fit2.wave, fit2.stki[x2,y2],fit2.stkq[x2,y2],fit2.stku[x2,y2],fit2.stkv[x1,y2]

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

	# 1st inversion
	label1 = '1st Inversion'
	if '-label1' in sys.argv:
		label1 = sys.argv[sys.argv.index("-label1")+1]
	# 2nd inversion
	label2 = '2nd Inversion'
	if '-label2' in sys.argv:
		label2 = sys.argv[sys.argv.index("-label2")+1]

	#############################################################
	#					PLOTTING STUFF					#
	#############################################################
	# Change to abs. wavelength to the line core of the first number
	if conf1['mode'] == "MC":
		ll0 = sir.determine_line_core(os.path.join(conf1["path"],conf1["line"]),num)

	label_x = 0

	temp = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if temp != '0':
		ll1 -= float(temp)
		ll2 -= float(temp)
		ll1_fit -= float(temp)
		ll2_fit -= float(temp)
		label_x = temp
	else:
		# Keep the same wavelengths as in the Grid file
		label_x = str(ll0)
		ll1 -= ll0
		ll2 -= ll0
		ll1_fit -= ll0
		ll2_fit -= ll0


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
		if "-one" not in sys.argv:
			llabel = "Syn. Profile ("
		else:
			llabel = "Syn. Profile"
	else:
		if "-one" not in sys.argv:
			llabel = "Obs. Profile ("
		else:
			llabel = "Obs. Profile"
	if "-one" in sys.argv:
		ax1.plot(ll1, I1, "x", label=llabel)
		ax2.plot(ll1, Q1, "x", label=llabel)
		ax3.plot(ll1, U1, "x", label=llabel)
		ax4.plot(ll1, V1, "x", label=llabel)
	else:
		ax1.plot(ll1, I1, "x", label=llabel + label1 + ")")
		ax1.plot(ll2, I2, "x", label=llabel + label2 + ")")
		ax2.plot(ll1, Q1, "x", label=llabel + label1 + ")")
		ax2.plot(ll2, Q2, "x", label=llabel + label2 + ")")
		ax3.plot(ll1, U1, "x", label=llabel + label1 + ")")
		ax3.plot(ll2, U2, "x", label=llabel + label2 + ")")
		ax4.plot(ll1, V1, "x", label=llabel + label1 + ")")
		ax4.plot(ll2, V2, "x", label=llabel + label2 + ")")

	ax1.plot(ll1, I_fit1, "-", label = "Best Fit (" + label1 + ")")
	ax1.plot(ll2, I_fit2, "-", label = "Best Fit (" + label2 + ")")
	ax2.plot(ll1, Q_fit1, "-", label = "Best Fit (" + label1 + ")")
	ax2.plot(ll2, Q_fit2, "-", label = "Best Fit (" + label2 + ")")
	ax3.plot(ll1, U_fit1, "-", label = "Best Fit (" + label1 + ")")
	ax3.plot(ll2, U_fit2, "-", label = "Best Fit (" + label2 + ")")
	ax4.plot(ll1, V_fit1, "-", label = "Best Fit (" + label1 + ")")
	ax4.plot(ll2, V_fit2, "-", label = "Best Fit (" + label2 + ")")

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
	
	
	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' [\AA]', loc='center')
	ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' [\AA]', loc='center')

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
		ax1.set_ylim(0.9*np.min(np.abs(np.append(I1,[I_fit1,I2,I_fit2]))), 1.1*np.max(np.abs(np.append(I1,[I_fit1,I2,I_fit2]))))
		ax2.set_ylim(-1.1*np.max(np.abs(np.append(Q1,[Q_fit1,Q2,Q_fit2]))), 1.1*np.max(np.abs(np.append(Q1,[Q_fit1,Q2,Q_fit2]))))
		ax3.set_ylim(-1.1*np.max(np.abs(np.append(U1,[U_fit1,U2,U_fit2]))), 1.1*np.max(np.abs(np.append(U1,[U_fit1,U2,U_fit2]))))
		ax4.set_ylim(-1.1*np.max(np.abs(np.append(V1,[V_fit1,V2,V_fit2]))), 1.1*np.max(np.abs(np.append(V1,[V_fit1,V2,V_fit2]))))
		ax1.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout(pad=2.5)

	plt.savefig(savepath + "inversion2_stokes_x1" + str(x1) + "_y1" + str(y1) +"_x2" + str(x2) + "_y2" + str(y2) + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-gamma", "-phi", "-z", "-Pg","-rho"]
	index  = [1,2,3,4,5,6,7,8,9,10]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$"]
	titles   = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength B", r"LOS Velocity $\mathrm{v}_{\mathrm{los}}$",
			 r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$"]

	i = 0

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
					llabel = "Best Fit M. 1 ("
				else:
					llabel = "Syn Model ("

				if "-one" not in sys.argv:
					ax1.plot(syn1.tau, syn1.get_attribute(inputs[i][1:])[x1,y1], label=f"{llabel}{label1})",color='#0C5DA5')
					ax1.plot(syn2.tau, syn2.get_attribute(inputs[i][1:])[x2,y2], label=f"{llabel}{label2})",color='#845B97')
				else:
					ax1.plot(syn1.tau, syn1.get_attribute(inputs[i][1:])[x1,y1], label=f"{llabel[:-2]})",color='#0C5DA5')
			if conf1['mode'] == "2C":
					llabel = "Best Fit M. 2 ("
			else:
					llabel = "Best Fit ("
			
			ax1.plot(fit1.tau, fit1.get_attribute(inputs[i][1:])[x1,y1], label=f"{llabel}{label1})",color = '#FF2C00')
			ax1.plot(fit2.tau, fit2.get_attribute(inputs[i][1:])[x2,y2], label=f"{llabel}{label2})",color='#00B945')

			if "-err" in sys.argv:
				if conf1['mode'] == "2C":
					# Error of fit
					ax1.fill_between(syn1.tau, syn1.get_attribute(inputs[i][1:])[x1,y1] - err11.get_attribute(inputs[i][1:])[x1,y1],
								syn1.get_attribute(inputs[i][1:])[x1,y1] + err11.get_attribute(inputs[i][1:])[x1,y1], alpha = 0.5,
								color='#0C5DA5', lw=0)
					ax1.fill_between(syn2.tau, syn2.get_attribute(inputs[i][1:])[x2,y2] - err22.get_attribute(inputs[i][1:])[x2,y2],
								syn2.get_attribute(inputs[i][1:])[x2,y2] + err22.get_attribute(inputs[i][1:])[x2,y2], alpha = 0.5,
								color='#845B97', lw=0)
				
				# Error of fit
				ax1.fill_between(fit1.tau, fit1.get_attribute(inputs[i][1:])[x1,y1] - err1.get_attribute(inputs[i][1:])[x1,y1],
							 fit1.get_attribute(inputs[i][1:])[x1,y1] + err1.get_attribute(inputs[i][1:])[x1,y1], alpha = 0.5,
							 color='#FF2C00', lw=0)
				ax1.fill_between(fit2.tau, fit2.get_attribute(inputs[i][1:])[x2,y2] - err2.get_attribute(inputs[i][1:])[x2,y2],
							 fit2.get_attribute(inputs[i][1:])[x2,y2] + err2.get_attribute(inputs[i][1:])[x2,y2], alpha = 0.5,
							 color='#00B945', lw=0)


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
			plt.savefig(savepath + "inversion2_x1" + str(x1) + "_y1"  + str(y1) + "_x2" + str(x2) + "_y2"  + str(y2) + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = syn1.tau[-1]
	if conf1["mode"] == "MC" or conf1['mode'] == "2C":
		syn1.set_limit(lim_max)
		syn1.set_limit(lim_max)
	phy1.set_limit(lim_max)
	err1.set_limit(lim_max)
	phy2.set_limit(lim_max)
	err2.set_limit(lim_max)

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	if conf1["mode"] == "MC" or conf1['mode'] == "2C":
		if conf1["mode"] == "2C":
			llabel = "Best Fit M. 1"
		else:
			llabel = "Syn. Model"
		if "-one" not in sys.argv:
			ax1.plot(syn1.tau, syn1.T[x1,y1], label=f"{llabel} ({label1})", color='#0C5DA5')
			ax2.plot(syn1.tau, syn1.B[x1,y1], label=f"{llabel} ({label1})", color='#0C5DA5')
			ax3.plot(syn1.tau, syn1.vlos[x1,y1], label=f"{llabel} ({label1})", color='#0C5DA5')
			ax4.plot(syn1.tau, syn1.gamma[x1,y1], label=f"{llabel} ({label1})", color='#0C5DA5')

		else:
			ax1.plot(syn1.tau, syn1.T[x1,y1], label=f"{llabel}", color='#0C5DA5')
			ax2.plot(syn1.tau, syn1.B[x1,y1], label=f"{llabel}", color='#0C5DA5')
			ax3.plot(syn1.tau, syn1.vlos[x1,y1], label=f"{llabel}", color='#0C5DA5')
			ax4.plot(syn1.tau, syn1.gamma[x1,y1], label=f"{llabel}", color='#0C5DA5')

	if conf1['mode'] == "2C":
		llabel = "Best Fit M. 2"
	else:
		llabel = "Best Fit"
	ax1.plot(phy1.tau, phy1.T[x1,y1], label=f"{llabel} ({label1})", color='#FF2C00')
	ax2.plot(phy1.tau, phy1.B[x1,y1], label=f"{llabel} ({label1})", color='#FF2C00')
	ax3.plot(phy1.tau, phy1.vlos[x1,y1], label=f"{llabel} ({label1})", color='#FF2C00')
	ax4.plot(phy1.tau, phy1.gamma[x1,y1], label=f"{llabel} ({label1})", color='#FF2C00')

	if conf1["mode"] == "MC" or conf1['mode'] == "2C":
		if conf1["mode"] == "2C":
			llabel = "Best Fit M. 2"
		else:
			llabel = "Syn. Model"
		if "-one" not in sys.argv:
			ax1.plot(syn2.tau, syn2.T[x2,y2], label=f"{llabel} ({label2})", color='#00B945')
			ax2.plot(syn2.tau, syn2.B[x2,y2], label=f"{llabel} ({label2})", color='#00B945')
			ax3.plot(syn2.tau, syn2.vlos[x2,y2], label=f"{llabel} ({label2})", color='#00B945')
			ax4.plot(syn2.tau, syn2.gamma[x2,y2], label=f"{llabel} ({label2})", color='#00B945')

	if conf1['mode'] == "2C":
		llabel = "Best Fit M. 2"
	else:
		llabel = "Best Fit"
	ax1.plot(phy2.tau, phy2.T[x2,y2], label=f"{llabel} ({label2})", color='#FF9500')
	ax2.plot(phy2.tau, phy2.B[x2,y2], label=f"{llabel} ({label2})", color='#FF9500')
	ax3.plot(phy2.tau, phy2.vlos[x2,y2], label=f"{llabel} ({label2})", color='#FF9500')
	ax4.plot(phy2.tau, phy2.gamma[x2,y2], label=f"{llabel} ({label2})", color='#FF9500')

	if "-err" in sys.argv:
		if conf1['mode'] == "2C":
			ax1.fill_between(syn1.tau, syn1.T[x1,y1] - err11.T[x1,y1],
					 syn1.T[x1,y1] + err11.T[x1,y1], alpha = 0.5,
					 color='#0C5DA5', lw=0)
			ax2.fill_between(syn1.tau, syn1.B[x1,y1] - err11.B[x1,y1],
						syn1.B[x1,y1] + err11.B[x1,y1], alpha = 0.5,
						color='#0C5DA5', lw=0)
			ax3.fill_between(syn1.tau, syn1.vlos[x1,y1] - err11.vlos[x1,y1],
						syn1.vlos[x1,y1] + err11.vlos[x1,y1], alpha = 0.5,
						color='#0C5DA5', lw=0)
			ax4.fill_between(syn1.tau, syn1.gamma[x1,y1] - err11.gamma[x1,y1],
						syn1.gamma[x1,y1] + err11.gamma[x1,y1], alpha = 0.5,
						color='#0C5DA5', lw=0)
			
			ax1.fill_between(syn2.tau, syn2.T[x2,y2] - err22.T[x2,y2],
					 syn2.T[x2,y2] + err22.T[x2,y2], alpha = 0.5,
					 color='#00B945', lw=0)
			ax2.fill_between(syn2.tau, syn2.B[x2,y2] - err22.B[x2,y2],
					syn2.B[x2,y2] + err22.B[x2,y2], alpha = 0.5,
					color='#00B945', lw=0)
			ax3.fill_between(syn2.tau, syn2.vlos[x2,y2] - err22.vlos[x2,y2],
					syn2.vlos[x2,y1] + err22.vlos[x2,y2], alpha = 0.5,
					color='#00B945', lw=0)
			ax4.fill_between(syn2.tau, syn2.gamma[x2,y2] - err22.gamma[x2,y2],
					syn2.gamma[x2,y2] + err22.gamma[x2,y2], alpha = 0.5,
					color='#00B945', lw=0)
			
	
	
		ax1.fill_between(phy1.tau, phy1.T[x1,y1] - err1.T[x1,y1],
					 phy1.T[x1,y1] + err1.T[x1,y1], alpha = 0.5,
					 color='#FF2C00', lw=0)
		ax2.fill_between(phy1.tau, phy1.B[x1,y1] - err1.B[x1,y1],
					 phy1.B[x1,y1] + err1.B[x1,y1], alpha = 0.5,
					 color='#FF2C00', lw=0)
		ax3.fill_between(phy1.tau, phy1.vlos[x1,y1] - err1.vlos[x1,y1],
					 phy1.vlos[x1,y1] + err1.vlos[x1,y1], alpha = 0.5,
					 color='#FF2C00', lw=0)
		ax4.fill_between(phy1.tau, phy1.gamma[x1,y1] - err1.gamma[x1,y1],
					 phy1.gamma[x1,y1] + err1.gamma[x1,y1], alpha = 0.5,
					 color='#FF2C00', lw=0)
		ax1.fill_between(phy2.tau, phy2.T[x2,y2] - err2.T[x2,y2],
					 phy2.T[x2,y2] + err2.T[x2,y2], alpha = 0.5,
					 color='#FF9500', lw=0)
		ax2.fill_between(phy2.tau, phy2.B[x2,y2] - err2.B[x2,y2],
					 phy2.B[x2,y2] + err2.B[x2,y2], alpha = 0.5,
					 color='#FF9500', lw=0)
		ax3.fill_between(phy2.tau, phy2.vlos[x2,y2] - err2.vlos[x2,y2],
					 phy2.vlos[x2,y2] + err2.vlos[x2,y2], alpha = 0.5,
					 color='#FF9500', lw=0)
		ax4.fill_between(phy2.tau, phy2.gamma[x2,y2] - err2.gamma[x2,y2],
					 phy2.gamma[x2,y2] + err2.gamma[x2,y2], alpha = 0.5,
					 color='#FF9500', lw=0)
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

	  
	plt.savefig(savepath + "inversion2_result_x1" + str(x1) + "_y1"  + str(y1) + "_x2" + str(x2) + "_y2"  + str(y2) + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf1 = sir.read_config(sys.argv[1], check=False)
	conf2 = sir.read_config(sys.argv[4], check=False)

	if conf1['path'][:2] == "./":
		conf1['path'] = sys.argv[1][:sys.argv[1].rfind('/')+1] + conf1['path'][2:]
	if conf2['path'][:2] == "./":
		conf2['path'] = sys.argv[4][:sys.argv[4].rfind('/')+1] + conf2['path'][2:]
	inversion_2(conf1, int(sys.argv[2]),int(sys.argv[3]),conf2, int(sys.argv[5]),int(sys.argv[6]))




