"""
Plots the result of the SIR inversion
"""

import numpy as np 
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(os.path.join(os.path.dirname(__file__), "../../tools"))
import sir
import model as m
import definitions as d
import matplotlib.pyplot as plt
from change_config_path import change_config_path
from os.path import exists
import matplotlib
import os


def help():
	"""
	Help Page
	"""
	print("plot_inversion - Plots the result of an inversion")
	print("Usage: python plot_inversion [OPTION]")
	print()
	sir.option("[1. Pos]","Config Model 1")
	sir.option("[2. Pos]","Number of the Model, Profile (starting at 1)")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
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
	sir.option("-line","Number of the line considered (Default: take first one)")
	sys.exit()


def inversion(conf, num):
	"""
	Plots the result of one inversion

	Parameter
	---------
	config : dict
		Config info
	num : int
		Model number
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
				
	# Check if path exists
	if not exists(conf['path']):
		Inp = input("[NOTE] Path does not exist. You want to overwrite it with the actual path? [y/n] ")
		if Inp == "y":
			change_config_path(conf, os.path.abspath(os.getcwd()))
	############################
	# READ INPUT AND LOAD DATA #
	############################
	path = conf["path"]
	# Which line is used (important for determination which npy file is loaded)
	if '-line' in sys.argv:
		line = int(sys.argv[sys.argv.index("-line")+1])
	else:
		Grid = sir.read_grid(os.path.join(conf['path'],d.Grid))
		line = int(Grid['Line'][0][0])

	syn = m.model(os.path.join(path, conf["model_out"]))		   				# Synthesis Models
	phy = m.model(os.path.join(path, conf["inv_out"] + d.inv_models))		# Inversion Models
	err = m.error(os.path.join(path, conf["inv_out"] + d.inv_errors))		# Inversion Models

	obs = np.load(os.path.join(path, conf["noise_out"] + f".npy"))[num-1]		# Synthesis Profiles
	fit = np.load(os.path.join(path, conf["inv_out"] + d.inv_stokes + ".npy"))[num-1]   # Inversion Profiles

	instrument = conf["instrument"]  # Instrument used for labels

	# Observation from synthesis
	ll, I, Q, U, V = obs[1], obs[2], obs[3], obs[4], obs[5]
	
	# Load fit data
	ll_fit, I_fit, Q_fit, U_fit, V_fit = fit[1],fit[2],fit[3],fit[4],fit[5s]

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = path + "/" + sys.argv[sys.argv.index("-save")+1]
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
	
	# Additional label
	add_label = '_'
	if '-label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]

	##################
	# PLOTTING STUFF #
	##################
	# Change to abs. wavelength to the line core of the first number
	Line = sir.read_line(conf["line"])
	Lines = Line['Line']
	ll0 = Line['wavelength'][np.where(Lines == line)[0][0]]
	ll /= 1e3  # from mA in A
	ll_fit /= 1e3  # from mA in A
	ll += ll0
	ll_fit += ll0

	label_x = 0
	# Change to  6300 as relative lambda0
	if instrument == "Hinode":
		ll -= 6300.0000
		ll_fit -= 6300.0000
		label_x = '6300'
	# Change to 15600 as relative lambda0
	elif instrument == "GRIS":
		ll -= 15600.0000
		ll_fit -= 15600.0000
		label_x = '15600'

	else:
		# Ask for the point to which it should be rel.
		temp = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
		if temp != '0':
			ll -= float(temp)
			ll_fit -= float(temp)
			label_x = temp
		else:
			# Keep the same wavelengths as in the Grid file
			label_x = str(ll0)
			ll -= ll0
			ll_fit -= ll0
	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]  # Get colors used in the actual cycle
	########################
	#  Plot I, Q, U and V  #
	########################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, 1, figsize=(12, 14), sharex=True, gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, figsize=(16, 12))

	# Add label
	ax2.scatter([], [], color="w", alpha=0, label=add_label)

	############################
	# Plot the Stokes profiles #
	############################
	llabel = "Synthesised Profile"
	ax1.plot(ll, I, "x", label=llabel)
	ax2.plot(ll, Q, "x", label=llabel)
	ax3.plot(ll, U, "x", label=llabel)
	ax4.plot(ll, V, "x", label=llabel)

	ax1.plot(ll, I_fit, "-", label="Best Fit")
	ax2.plot(ll, Q_fit, "-", label="Best Fit")
	ax3.plot(ll, U_fit, "-", label="Best Fit")
	ax4.plot(ll, V_fit, "-", label="Best Fit")

	ax1.plot(ll, I, "-", alpha=0.3, color=colors[0])
	ax2.plot(ll, Q, "-", alpha=0.3, color=colors[0])
	ax3.plot(ll, U, "-", alpha=0.3, color=colors[0])
	ax4.plot(ll, V, "-", alpha=0.3, color=colors[0])

	# Set x limits
	ax1.set_xlim(ll[0], ll[-1])
	ax2.set_xlim(ll[0], ll[-1])
	ax3.set_xlim(ll[0], ll[-1])
	ax4.set_xlim(ll[0], ll[-1])

	#######################################################################
	# Set labels												          #
	# The labels depend on the arguments and the chosen line		      #
	# The code is so long as I tried to make it as automatic as possible  #
	#######################################################################
	ax1.set_ylabel(r'$I / I_c$')
	ax2.set_ylabel(r'$Q / I_c$')
	ax3.set_ylabel(r'$U / I_c$')
	ax4.set_ylabel(r'$V / I_c$')

	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' $[\AA]$')
	ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' $[\AA]$')

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
		ax1.set_ylim(0.9*np.min(np.abs(np.append(I,I_fit))), 1.1*np.max(np.abs(np.append(I,I_fit))))
		ax2.set_ylim(-1.1*np.max(np.abs(np.append(Q,Q_fit))), 1.1*np.max(np.abs(np.append(Q,Q_fit))))
		ax3.set_ylim(-1.1*np.max(np.abs(np.append(U,U_fit))), 1.1*np.max(np.abs(np.append(U,U_fit))))
		ax4.set_ylim(-1.1*np.max(np.abs(np.append(V,V_fit))), 1.1*np.max(np.abs(np.append(V,V_fit))))
		ax1.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout(pad=2.5)

	plt.savefig(savepath + "inversion_stokes" + str(num) + add)

	############################
	# Plot physical parameters #
	############################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho"]
	index = [1,2,3,4,5,6,7,8,9,10]
	att = ["T", "Pe", "vmicro", "B", "vlos", "gamma", "phi", "z", "Pg", "rho"] # For getting the value from the class
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$",
			  r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]",
			  r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]",
			  r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$",
			  r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$"
			  ]
	titles = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$",
			  r"Magnetic Field Strength B", r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
			  r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$"]
	err_limits = [1e4, 1e3, 10, 1e4, 1e5, 360, 360 ,1e3, 1e6, 1e-3]	# Limits were the error bar is not plotted to really see the results

	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Get colors used in the actual cycle

	# Limit the plots to log_tau -3
	lim_max = -3
	syn.set_limit(lim_max)
	phy.set_limit(lim_max)
	err.set_limit(lim_max)

	# Correct phi range
	syn.correct_phi()
	phy.correct_phi()

	if syn.log_tau[num-1,-1] > lim_max:
		lim_max = syn.log_tau[num-1,-1]

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8,7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)
			# Plot
			# Synthesised model = Real model
			ax1.plot(syn.log_tau[num-1], syn.get_attribute(att[i])[num-1], label="Model")

			ax1.plot(phy.log_tau[num-1], phy.get_attribute(att[i])[num-1], label="Best Fit Model", color=colors[1])

			if err.get_attribute(att[i])[num-1,0] > err_limits[i]:
				print(f'Note that the error for {att[i]} is very big with a value of {err.get_attribute(att[i])[num-1,0]} and is not plotted.')
			else:
				# Error of fit
				ax1.fill_between(phy.log_tau[num-1], phy.get_attribute(att[i])[num-1] - err.get_attribute(att[i])[num-1],
								 phy.get_attribute(att[i])[num-1] + err.get_attribute(att[i])[num-1], alpha = 0.5,
								 color=colors[1], lw=0)

			# Set x limits
			ax1.set_xlim(phy.log_tau[num-1,0], phy.log_tau[num-1,-1])

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
			plt.savefig(savepath + "inversion_" + str(inputs[i][1:]) + str(num) + add)
		
	# Plot T,B,vlos, inc in one figure

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	llabel = "Synthesised Model"
	ax1.plot(syn.log_tau[num-1], syn.T[num-1], label=llabel, color='#0C5DA5')
	ax2.plot(syn.log_tau[num-1], syn.B[num-1], label=llabel, color='#00B945')
	ax3.plot(syn.log_tau[num-1], syn.vlos[num-1], label=llabel, color='#845B97')
	ax4.plot(syn.log_tau[num-1], syn.gamma[num-1], label=llabel, color='#054907')

	llabel = "Best Fit"
	ax1.plot(phy.log_tau[num-1], phy.T[num-1], label=llabel, color='#FF2C00')
	ax2.plot(phy.log_tau[num-1], phy.B[num-1], label=llabel, color='#FF9500')
	ax3.plot(phy.log_tau[num-1], phy.vlos[num-1], label=llabel, color='#474747')
	ax4.plot(phy.log_tau[num-1], phy.gamma[num-1], label=llabel, color='#C20078')

	if err.T[num-1,0] > err_limits[0]:
		print(f'Note that the error for T is very big with a value of {err.T[num-1,0]} K and is not plotted.')
	else:
		ax1.fill_between(phy.log_tau[num-1], phy.T[num-1] - err.T[num-1],
				 phy.T[num-1] + err.T[num-1], alpha = 0.5,
				 color='#FF2C00', lw=0)
	if err.B[num-1,0] > err_limits[3]:
		print(f'Note that the error for B is very big with a value of {err.B[num-1,0]} G and is not plotted.')
	else:
		ax2.fill_between(phy.log_tau[num-1], phy.B[num-1] - err.B[num-1],
				 phy.B[num-1] + err.B[num-1], alpha = 0.5,
				 color='#FF9500', lw=0)
	if err.vlos[num-1,0] > err_limits[4]:
		print(f'Note that the error for vlos is very big with a value of {err.vlos[num-1,0]} km/s and is not plotted.')
	else:
		ax3.fill_between(phy.log_tau[num-1], phy.vlos[num-1] - err.vlos[num-1],
				 phy.vlos[num-1] + err.vlos[num-1], alpha = 0.5,
				 color='#474747', lw=0)
	if err.gamma[num-1,0] > err_limits[5]:
		print(f'Note that the error for gamma is very big with a value of {err.gamma[num-1,0]} and is not plotted.')
	else:
		ax4.fill_between(phy.log_tau[num-1], phy.gamma[num-1] - err.gamma[num-1],
				 phy.gamma[num-1] + err.gamma[num-1], alpha = 0.5,
				 color='#C20078', lw=0)

	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(phy.log_tau[num-1,0], lim_max)
	ax2.set_xlim(phy.log_tau[num-1,0], lim_max)
	ax3.set_xlim(phy.log_tau[num-1,0], lim_max)
	ax4.set_xlim(phy.log_tau[num-1,0], lim_max)

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
		ax3.set_xlabel(r"$\log \tau_{500}$")
	ax4.set_xlabel(r"$\log \tau_{500}$")

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

	  
	plt.savefig(savepath + "inversion_result" + str(num) + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])
	inversion(conf, int(sys.argv[2]))

