"""
Plots the result of the SIR synthesis
"""

import numpy as np 
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../..")) 
import sir, model as m
import definitions as d
import matplotlib.pyplot as plt
from os.path import exists
import os

# TODO implement the class model/error

def help():
	"""
	Help Page
	"""
	print("plot_inversion_2 - Plots the result of two inversions")
	print("Usage: python plot_inversion_2 [OPTION]")
	print()
	sir.option("[1. Pos]","Config Model 1")
	sir.option("[2. Pos]","Number of the Model, Profile 1")
	sir.option("[3. Pos]","Config Model 2")
	sir.option("[4. Pos]","Number of the Model, Profile 2")

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
	sir.option("-inc","Plot inclination by subtracting in deg")
	sir.option("-azi","Plot azimuth by adding in deg")
	sir.option("-z","Plot real height in km")
	sir.option("-Pg","Plot gas pressure in dyn/cm^2")
	sir.option("-rho","Plot density")
	sir.option("-syn","Synthesised model .mod file")
	sir.option("-vertical","Plot spectra vertically")
	sir.option("-line:","Number of the line considered (Default: take first one)")
	sys.exit()

def inversion_2(conf1, num1, conf2, num2):
	""""
	Plots the result of one inversion

	Parameter
	---------
	conf1 : dict
		Config infos for 1st inversion
	num1 : int
		Model number for 1st inversion
	conf2 : dict
		Config infos for 2nd inversion
	num2 : int
		Model number for 2nd inversion
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

	
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################

	path1 = conf1["path"]
	syn1 = np.load(os.path.join(path1,conf1["model_out"]))[num1]				# Synthesis Models 1
	phy1 = np.load(os.path.join(path1,conf1["inv_out"] + d.end_models))[num1]		# Inversion Models 1
	err1 = np.load(os.path.join(path1,conf1["inv_out"] + d.end_errors))[num1]		# Inversion Models 1
	obs1 = np.load(os.path.join(path1,conf1["syn_out"]))[num1]	# Synthesis Profiles 1
	fit1 = np.load(os.path.join(path1,conf1["inv_out"] + d.end_stokes))[num1]	# Inversion Profiles 1
	path2 = conf2["path"]
	syn2 = np.load(os.path.join(path2,conf2["model_out"]))[num2]				# Synthesis Models 2
	phy2 = np.load(os.path.join(path2,conf2["inv_out"] + d.end_models))[num2]		# Inversion Models 2
	err2 = np.load(os.path.join(path2,conf2["inv_out"] + d.end_errors))[num2]		# Inversion Models 2
	obs2 = np.load(os.path.join(path2,conf2["syn_out"]))[num2]	# Synthesis Profiles 2
	fit2 = np.load(os.path.join(path2,conf2["inv_out"] + d.end_stokes))[num2]	# Inversion Profiles 2

	
	instrument1 = conf1["instrument"] # Instrument used for labels
	instrument2 = conf1["instrument"] # Instrument used for labels

	if instrument1 != instrument2:
		print("[WARN] Script should only be used for the same instrument!")
		sys.exit()
	else:
		instrument = instrument1

	# Observation from synthesis
	ll1, I1, Q1, U1, V1 = obs1[1],obs1[2],obs1[3],obs1[4],obs1[5]
	ll2, I2, Q2, U2, V2 = obs2[1],obs2[2],obs2[3],obs2[4],obs2[5]
	
	# Load fit data
	ll_fit1, I_fit1, Q_fit1, U_fit1, V_fit1 = fit1[1],fit1[2],fit1[3],fit1[4],fit1[5]
	ll_fit2, I_fit2, Q_fit2, U_fit2, V_fit2 = fit2[1],fit2[2],fit2[3],fit2[4],fit2[5]
	
	# vlos in km /s
	syn1[:,:,5] = syn1[:,:,5]/1e5
	phy1[:,:,5] = phy1[:,:,5]/1e5
	err1[:,:,5] = err1[:,:,5]/1e5
	syn2[:,:,5] = syn2[:,:,5]/1e5
	phy2[:,:,5] = phy2[:,:,5]/1e5
	err2[:,:,5] = err2[:,:,5]/1e5

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
	label1 = '1st Model'
	if '-label1' in sys.argv:
		label1 = sys.argv[sys.argv.index("-label1")+1]
	# 2nd inversion
	label2 = '2nd Model'
	if '-label2' in sys.argv:
		label2 = sys.argv[sys.argv.index("-label2")+1]

	#############################################################
	#					PLOTTING STUFF					#
	#############################################################
	# Change to abs. wavelength to the line core of the first number
	Line = sir.read_line(conf1["line"])
	Lines = Line['Line']
	ll0 = Line['wavelength'][np.where(Lines == line)[0][0]]
	ll += ll0
	ll_fit += ll0

	label_x = 0
	# Change to  6300 as relative lambda0
	if instrument == "Hinode" or instrument == "GRIS":
		ll1 -= d.ll_relative[instrument]
		ll2 -= d.ll_relative[instrument]
		ll1_fit -= d.ll_relative[instrument]
		ll2_fit -= d.ll_relative[instrument]
		label_x = str(d.ll_relative[instrument])
	
	else:
		# Ask for the point to which it should be rel.
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
	llabel = "Syn. Profile ("
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
		if instrument == "GRIS":
			fig.suptitle("Near-Infrared Lines", y=0.98, x=xtitle1)
		elif instrument == "Hinode":
			fig.suptitle("Visible Lines", y=0.98, x=xtitle1)
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

	plt.savefig(savepath + "inversion2_stokes" + n + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho"]
	index  = [1,2,3,4,5,6,7,8,9,10]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$"]
	titles   = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength B", r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
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
			ax1.plot(syn1[0,0], syn1[index[i],0], label=f"Syn. Model ({label1})",color='#0C5DA5')
			ax1.plot(syn2[0,0], syn2[index[i],0], label=f"Syn. Model ({label2})",color='#845B97')

			ax1.plot(fit1[0,0], fit1[index[i],0], label=f"Best Fit ({label1})",color = '#FF2C00')
			ax1.plot(fit2[0,0], fit2[index[i],0], label=f"Best Fit ({label2})",color='#00B945')

			# Error of fit
			ax1.fill_between(fit1[0,0], fit1[index[i],0] - err1[index[i],0],
						 fit1[index[i],0] + err1[index[i],0], alpha = 0.5,
						 color='#FF2C00', lw=0)
			ax1.fill_between(fit2[0,0], fit2[index[i],0] - err2[index[i],0],
						 fit2[index[i]] + err2[index[i]], alpha = 0.5,
						 color='00B945', lw=0)


			# Set xlimits
			ax1.set_xlim(fit1[00,][0], fit1[0][-1])

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{500}$")
			ax1.set_ylabel(labels[index[i]])

			if index[i] == 2 or index[i] == 9:
				ax1.semilogy()

			# Legend
			ax1.legend()
		
			ax1.set_title(titles[i])
			# set the spacing between subplots
			plt.tight_layout(pad=2)
			plt.savefig(savepath + "inversion2_" + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = syn1[0,0,-1]
	if lim_max < -3:
		lim_max = -3
		# Cut data so that the plot limits are adjusted to the shorted range
		I = np.where(syn1[0] < -3)[0][0]
		syn1 = syn1[:,0,0:I]
		syn2 = syn2[:,0,0:I]
		phy1 = phy1[:,0,0:I]
		err1 = err1[:,0,0:I]
		phy2 = phy2[:,0,0:I]
		err2 = err2[:,0,0:I]

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	llabel = "Syn. Model"
	ax1.plot(syn1[0,0], syn1[1,0], label=f"{llabel} ({label1})", color='#0C5DA5')
	ax2.plot(syn1[0,0], syn1[4,0], label=f"{llabel} ({label1})", color='#0C5DA5')
	ax3.plot(syn1[0,0], syn1[5,0], label=f"{llabel} ({label1})", color='#0C5DA5')
	ax4.plot(syn1[0,0], syn1[6,0], label=f"{llabel} ({label1})", color='#0C5DA5')

	llabel = "Best Fit"
	ax1.plot(syn1[0,0], phy1[1,0], label=f"{llabel} ({label1})", color='#FF2C00')
	ax2.plot(syn1[0,0], phy1[4,0], label=f"{llabel} ({label1})", color='#FF2C00')
	ax3.plot(syn1[0,0], phy1[5,0], label=f"{llabel} ({label1})", color='#FF2C00')
	ax4.plot(syn1[0,0], phy1[6,0], label=f"{llabel} ({label1})", color='#FF2C00')

	label = "Syn. Model"
	ax1.plot(syn2[0,0], syn2[1,0], label=f"{llabel} ({label2})", color='#00B945')
	ax2.plot(syn2[0,0], syn2[4,0], label=f"{llabel} ({label2})", color='#00B945')
	ax3.plot(syn2[0,0], syn2[5,0], label=f"{llabel} ({label2})", color='#00B945')
	ax4.plot(syn2[0,0], syn2[6,0], label=f"{llabel} ({label2})", color='#00B945')

	llabel = "Best Fit"
	ax1.plot(syn2[0,0], phy2[1,0], label=f"{llabel} ({label2})", color='#FF9500')
	ax2.plot(syn2[0,0], phy2[4,0], label=f"{llabel} ({label2})", color='#FF9500')
	ax3.plot(syn2[0,0], phy2[5,0], label=f"{llabel} ({label2})", color='#FF9500')
	ax4.plot(syn2[0,0], phy2[6,0], label=f"{llabel} ({label2})", color='#FF9500')

	ax1.fill_between(phy1[0,0], phy1[1,0] - err1[1,0],
				 phy1[1,0] + err1[1,0], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax2.fill_between(phy1[0,0], phy1[4,0] - err1[4,0],
				 phy1[4,0] + err1[4,0], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax3.fill_between(phy1[0,0], phy1[5,0] - err1[5,0],
				 phy1[5,0] + err1[5,0], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax4.fill_between(phy1[0,0], phy1[6,0] - err1[6,0],
				 phy1[6,0] + err1[6,0], alpha = 0.5,
				 color='#FF2C00', lw=0)
	ax1.fill_between(phy2[0,0], phy2[1,0] - err2[1,0],
				 phy2[1,0] + err2[1,0], alpha = 0.5,
				 color='#FF9500', lw=0)
	ax2.fill_between(phy2[0,0], phy2[4,0] - err2[4,0],
				 phy2[4,0] + err2[4,0], alpha = 0.5,
				 color='#FF9500', lw=0)
	ax3.fill_between(phy2[0,0], phy2[5,0] - err2[5,0],
				 phy2[5,0] + err2[5,0], alpha = 0.5,
				 color='#FF9500', lw=0)
	ax4.fill_between(phy2[0,0], phy2[6,0] - err2[6,0],
				 phy2[6,0] + err2[6,0], alpha = 0.5,
				 color='#FF9500', lw=0)
	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(syn1[0,0,0], lim_max)
	ax2.set_xlim(syn1[0,0,0], lim_max)
	ax3.set_xlim(syn1[0,0,0], lim_max)
	ax4.set_xlim(syn1[0,0,0], lim_max)

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
		if instrument == "GRIS":
			fig.suptitle("Near-Infrared Lines", y=0.98, x=xtitle2)
		elif instrument == "Hinode":
			fig.suptitle("Visible Lines", y=0.98, x=xtitle2)
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle2)

	if "-vertical" in sys.argv:	
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		# set the spacing between subplots
		plt.tight_layout(pad=2)

	  
	plt.savefig(savepath + "inversion2_result_" + str(num1) + "_"  + str(num2)+ add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf1 = sir.read_config(sys.argv[1])
	conf2 = sir.read_config(sys.argv[3])
	inversion_2(conf1, int(sys.argv[2]),conf2, int(sys.argv[4]))




