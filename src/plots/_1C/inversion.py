"""
Plots the result of the SIR synthesis
"""

import numpy as np 
import sys, os
sys.path.append(sys.path[0] + "/../..")
sys.path.append(sys.path[0] + "/../tools")
import sir, obs
import definitions as d
import matplotlib.pyplot as plt
from os.path import exists
import matplotlib
import os
from change_config_path import change_config_path

def help():
	"""
	Help Page
	"""
	print("plot_inversion - Plots the result of an inversion")
	print("Usage: python plot_inversion [OPTION]")
	print()
	sir.option("[1. Pos.]","Config File")
	sir.option("[2. Pos.]","x position (of the complete map)")
	sir.option("[3. Pos.]","y position (of the complete map)")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
	sir.option("-title","Title in plots")
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
	sir.option("-xlim","Limits for the Stokes plot in A as a list")
	
	sys.exit()

def inversion(conf, x, y):
	""""
	Plots the result of one inversion

	Parameter
	---------
	config : dict
		Config infos
	x : int
		x position
	y : int
		y position
	"""
	# Import matplotlib library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/../mml.mplstyle'):
		plt.style.use(dirname + '/../mml.mplstyle')
	elif "mml" in plt.style.available:
		plt.style.use('mml')
	# Check if path exists
	if not exists(conf['path']):
		Inp = input("[NOTE] Path does not exist. You want to overwrite it with the actual path? [y/n] ")
		if Inp == "y":
			change_config_path(conf,os.path.abspath(os.getcwd()))
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################

	# Which line is used (important for determination which npy file is loaded)
	xlim = (None,None)
	if '-xlim' in sys.argv:
		xlim = [float(i) for i in sys.argv[sys.argv.index("-xlim")+1].split(',')]

	path = conf["path"]
	Map = conf['map']
	waves = np.load(os.path.join(path, conf['waves']))
	waves_inv = np.load(os.path.join(path,conf['inv_out']) + d.end_wave)
	range_wave1 = sir.angstrom_to_pixel(waves, conf['range_wave'])
	range_wave2 = sir.angstrom_to_pixel(waves_inv, conf['range_wave'])
	phy = np.load(os.path.join(path, conf["inv_out"] + d.end_models))[x - Map[0],y - Map[2]]    # Inversion Models
	err = np.load(os.path.join(path, conf["inv_out"] + d.end_errors))[x - Map[0],y - Map[2]]    # Inversion Errors
	obs1 = obs.load_data(conf, filename=conf['cube_inv'])[x, y]  # Observation Profiles
	fit = np.load(os.path.join(path, conf["inv_out"] + d.end_stokes))[x - Map[0],y - Map[2]]   # Inversion Profiles

	instrument = conf["instrument"] # Instrument used for labels

	# Observation from synthesis
	I, Q, U, V = obs1[0], obs1[1], obs1[2], obs1[3]
	
	# Load fit data
	I_fit, Q_fit, U_fit, V_fit = fit[0], fit[1], fit[2], fit[3]
	ll = np.copy(waves)

	# vlos in km /s
	phy[5] = phy[5]/1e5
	err[5] = err[5]/1e5

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = path + "/" + sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath[:savepath.rfind('/')]):
			os.mkdir(savepath[:savepath.rfind('/')])
	
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

	#############################################################
	# PLOTTING STUFF #
	#############################################################

	label_x = 0
	# Change to relative wavelength if wished ( see definitions.py file)
	if instrument in d.ll_relative:
		ll -= d.ll_relative[instrument]
		if xlim != (None, None):
			xlim[0] -= d.ll_relative[instrument]
			xlim[1] -= d.ll_relative[instrument]
		label_x = str(d.ll_relative[instrument])

	else:
		# Ask for the point to which it should be rel.
		temp = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
		if temp != '0':
			ll -= float(temp)
			label_x = temp

	# Delete not fitted part
	ll_temp = np.copy(ll)
	I_temp = np.copy(I_fit)
	Q_temp = np.copy(Q_fit)
	U_temp = np.copy(U_fit)
	V_temp = np.copy(V_fit)
	ll_fit = []
	I_fit = []
	Q_fit = []
	U_fit = []
	V_fit = []
	for i in range(len(range_wave2)):
		ll_fit.append(ll_temp[range_wave2[i][0]:range_wave2[i][1]+1])
		I_fit.append(I_temp[range_wave2[i][0]:range_wave2[i][1]+1])
		Q_fit.append(Q_temp[range_wave2[i][0]:range_wave2[i][1]+1])
		U_fit.append(U_temp[range_wave2[i][0]:range_wave2[i][1]+1])
		V_fit.append(V_temp[range_wave2[i][0]:range_wave2[i][1]+1])

	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Get colors used in the actual cycle

	########################
	#  Plot I, Q, U and V  #
	########################
	if "-vertical" in sys.argv:
		fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 14), sharex=True,
													gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

	# Add label
	ax2.scatter([], [], color="w", alpha=0, label=add_label)

	############################
	# Plot the Stokes profiles #
	############################
	llabel = "Observation"
	ax1.plot(ll, I, ".", label=llabel)
	ax2.plot(ll, Q, ".", label=llabel)
	ax3.plot(ll, U, ".", label=llabel)
	ax4.plot(ll, V, ".", label=llabel)
	ax1.plot(ll, I, "-", alpha=0.5, color=colors[0])
	ax2.plot(ll, Q, "-", alpha=0.5, color=colors[0])
	ax3.plot(ll, U, "-", alpha=0.5, color=colors[0])
	ax4.plot(ll, V, "-", alpha=0.5, color=colors[0])

	for i in range(len(ll_fit)):
		if i > 0:
			label1 = '_'
		else:
			label1= "Best Fit"
		ax1.plot(ll_fit[i], I_fit[i], "-", label=label1, color=colors[1])
		ax2.plot(ll_fit[i], Q_fit[i], "-", label=label1, color=colors[1])
		ax3.plot(ll_fit[i], U_fit[i], "-", label=label1, color=colors[1])
		ax4.plot(ll_fit[i], V_fit[i], "-", label=label1, color=colors[1])

	# Set xlimits
	ax1.set_xlim(xlim)
	ax2.set_xlim(xlim)
	ax3.set_xlim(xlim)
	ax4.set_xlim(xlim)

	#########################################################################
	# Set labels															#
	# The labels depend on the arguments and the chosen line				#
	# The code is so long as I tried to make it as automatic as possible	#
	#########################################################################
	ax1.set_ylabel(r'$I / I_c$')
	ax2.set_ylabel(r'$Q / I_c$')
	ax3.set_ylabel(r'$U / I_c$')
	ax4.set_ylabel(r'$V / I_c$')

	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' [\AA]')
	ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' [\AA]')

	#####################################################################
	# Set title															#
	# The position is relative to the chosen plot (vertical or not)  	#
	#####################################################################

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
		# ax1.set_ylim(0.9*np.min(np.abs(np.append(I,I_fit))), 1.1*np.max(np.abs(np.append(I,I_fit))))
		# ax2.set_ylim(-1.1*np.max(np.abs(np.append(Q,Q_fit))), 1.1*np.max(np.abs(np.append(Q,Q_fit))))
		# ax3.set_ylim(-1.1*np.max(np.abs(np.append(U,U_fit))), 1.1*np.max(np.abs(np.append(U,U_fit))))
		# ax4.set_ylim(-1.1*np.max(np.abs(np.append(V,V_fit))), 1.1*np.max(np.abs(np.append(V,V_fit))))
		ax1.legend(bbox_to_anchor=(1.01, 0.95))
		# set the spacing between subplots
		plt.tight_layout(pad=2, h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout(pad=2.5)

	plt.savefig(savepath + "inversion_stokes" + f"_{x}_{y}" + add)

	#############################
	# Plot physical parameters	#
	#############################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg", "-rho"]
	index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$",
				r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]",
				r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]",
				r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$",
				r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$"]
	titles = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$",
				r"Magnetic Field Strength B", r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
				r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$"]

	i = 0

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8, 7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)

			# Plot
			ax1.plot(phy[0], phy[index[i]], label="Best Fit Model", color='#FF2C00')

			# Error of fit
			ax1.fill_between(phy[0], phy[index[i]] - err[index[i]], phy[index[i]] + err[index[i]], alpha=0.5,
								color='#FF2C00', lw=0)

			# Set xlimits
			ax1.set_xlim(phy[0][0], phy[0][-1])

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
			plt.savefig(savepath + "inversion_" + str(inputs[i][1:]) + f"_{x}_{y}" + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = phy[0, -1]
	if lim_max < -3:
		lim_max = -3
		# Cut data so that the plot limits are adjusted to the shorted range
		n = np.where(phy[0] < -3)[0][0]
		err = err[:, 0:n]
		phy = phy[:, 0:n]
	if "-vertical" in sys.argv:
		fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(12, 16), sharex=True,
													gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12), sharex=True)

	llabel = "Best Fit"
	ax1.plot(phy[0], phy[1], label=llabel, color='#FF2C00')
	ax2.plot(phy[0], phy[4], label=llabel, color='#FF9500')
	ax3.plot(phy[0], phy[5], label=llabel, color='#474747')
	ax4.plot(phy[0], phy[6], label=llabel, color='#C20078')

	ax1.fill_between(phy[0], phy[1] - err[1], phy[1] + err[1], alpha=0.5, color='#FF2C00', lw=0)
	ax2.fill_between(phy[0], phy[4] - err[4], phy[4] + err[4], alpha=0.5, color='#FF9500', lw=0)
	ax3.fill_between(phy[0], phy[5] - err[5], phy[5] + err[5], alpha=0.5, color='#474747', lw=0)
	ax4.fill_between(phy[0], phy[6] - err[6], phy[6] + err[6], alpha=0.5, color='#C20078', lw=0)
	#################
	# Set limits	#
	#################
	ax1.set_xlim(phy[0][0], lim_max)
	ax2.set_xlim(phy[0][0], lim_max)
	ax3.set_xlim(phy[0][0], lim_max)
	ax4.set_xlim(phy[0][0], lim_max)

	Min, Max = ax2.get_ylim()
	ax2.set_ylim(Min, Max*1.15)
	Min, Max = ax3.get_ylim()
	ax3.set_ylim(Min, Max*1.15)
	Min, Max = ax4.get_ylim()
	ax4.set_ylim(Min, Max*1.15)

	#################
	# Set labels	#
	#################
	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r"$\log \tau_{500}$")
	ax4.set_xlabel(r"$\log \tau_{500}$")

	ax1.set_ylabel(labels[1])
	ax2.set_ylabel(labels[4])
	ax3.set_ylabel(labels[5])
	ax4.set_ylabel(labels[6])

	########################
	# Set title and legend #
	########################
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
		plt.tight_layout(pad=2, h_pad=0.0)
	else:
		# set the spacing between subplots
		plt.tight_layout(pad=2)

	plt.savefig(savepath + "inversion_result" + f"_{x}_{y}" + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])
	inversion(conf, int(sys.argv[2]), int(sys.argv[3]))

