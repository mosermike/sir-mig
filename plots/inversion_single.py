"""
Plots the result of the SIR synthesis
"""

import numpy as np 
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../src")) 
import sir
import definitions as d
import matplotlib.pyplot as plt
from os.path import exists

def _help():
	'''
	Help Page
	'''
	print("plot_inversion_single - Plots the result of a inversion without using a config file")
	print("Usage: python plot_inversion_single [OPTION]")
	print()
	sir.option("[1. Pos]","Mode (`1C`, `2C` or `MC`)")
	sir.option("[2. Pos]","Best Fit Profile")
	sir.option("[3. Pos]","Synthesis Profile")
	sir.option("[4. Pos]","Best Fit Model")
	sir.option("[5. Pos]","Synthesis Model")
	print()
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
	sir.option("-num","Number of the line considered")
	sir.option("-line","Line file if -num is used")
	sir.option("-instrument [str]","which instrument is used for relative wavelength (Hinode or GRIS at 15648A) (Mode `1C` and `2C`), otherwise asked")
	sys.exit()


def inversion_single_1C(fit : str, obs : str, phy : str):
	r"""
	
	Parameters
	----------
	fit : str
		Best Fit Profile
	obs : str
		Observation profiles
	phy : str
		Best Fit Model

	Returns
	-------
	None

	Other Parameters
	----------------
	Additional parameters given as an argument when the script is executed.

	-save [str]
		Additional save path (optional, default './')
	-add [str]
			Additional text in filenames (optional)
	-label [str]
		Add label text (optional)
	-title [str]
		Title in Stokes plot
	-T
		Plot temperature in K
	-Pe
		Plot electron pressure in dyn/cm^2
	-vmicro
		Plot microturbulence in cm/s
	-B
		Plot magnetic field strength in Gauss
	-vlos
		Plot line of sight velocity in cm/s
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
	-syn
		Synthesised model .mod file
	-vertical
		Plot spectra vertically
	-num [int]
		Number of the line considered
	-line [str]
		Line file if -num is used
	-instrument [str]
		which instrument is used (GRIS 15648A line or Hinode) for relative wavelength, otherwise asked
	
	"""

	# Import library
	sir.mpl_library()
	

	num = 0
	if '-num' in sys.argv:
		num = int(sys.argv[sys.argv.index("-num")+1])
	if '-line' in sys.argv:
		line = sys.argv[sys.argv.index("-line")+1]
	instrument =''
	if '-instrument' in sys.argv:
		instrument = sys.argv[sys.argv.index("-instrument")+1]

	# Observation from synthesis	
	ll, I, Q, U, V = sir.read_profile(obs, num)

	# Load fit data
	ll_fit, I_fit, Q_fit, U_fit, V_fit = sir.read_profile(fit, num)

	# Load physical parameters
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(phy)
	pars = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])

	# Load error of best fit model
	err = fit.replace(".per",".err")
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(err)
	pars_err = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])


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

	# Additional label
	add_label = '_'
	if '-label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]

	#############################################################
	#					PLOTTING STUFF					#
	#############################################################

	# Change rel. to 6300 for Hinode lines
	if '-num' in sys.argv:
		Line = sir.read_line(line)
		Lines = Line['Line']
		ll0 = Line['wavelength']
		ll0 = ll0[np.where(Lines == num)[0][0]]

	elif instrument == "Hinode":
		ll -= -6301.5012+6300.0000
		ll_fit -= -6301.5012+6300.0000
	elif instrument == "GRIS":
		ll -= -15648.514+15600.0000
		ll_fit -= -15648.514+15600.0000

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
	llabel = "Observation"
	ax1.plot(ll, I, "x", label=llabel)
	ax2.plot(ll, Q, "x", label=llabel)
	ax3.plot(ll, U, "x", label=llabel)
	ax4.plot(ll, V, "x", label=llabel)

	ax1.plot(ll, I_fit, "-", label = "Best Fit")
	ax2.plot(ll, Q_fit, "-", label = "Best Fit")
	ax3.plot(ll, U_fit, "-", label = "Best Fit")
	ax4.plot(ll, V_fit, "-", label = "Best Fit")

	# Set xlimits
	ax1.set_xlim(ll[0], ll[-1])
	ax2.set_xlim(ll[0], ll[-1])
	ax3.set_xlim(ll[0], ll[-1])
	ax4.set_xlim(ll[0], ll[-1])

	#######################################################################
	# Set labels												#
	# The labels depend on the arguments and the chosen line			#
	# The code is so long as I tried to make it as automatic as possible	#
	#######################################################################
	ax1.set_ylabel(r'$I / I_c$')
	ax2.set_ylabel(r'$Q / I_c$')
	ax3.set_ylabel(r'$U / I_c$')
	ax4.set_ylabel(r'$V / I_c$')
	# If -line and -num in arguments => determine rel. wavelength from there
	# Do not include if hinode lines are selected because then the rel. is 6300
	# to make it easier to be read from the plot
	title1 = None
	if '-num' in sys.argv:
		Line = sir.read_line(line)
		Lines = Line['Line']
		ll0 = Line['wavelength']
		ll0 = ll0[np.where(Lines == num)[0][0]]
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + str(ll0) + r' $[\AA]$', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - $' + str(ll0) + r' $[\AA]$', loc='center')
		title1 = "Fe I " + str(ll0)

	# Use the four 1.5 lines
	elif instrument == "GRIS":
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - 15600$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - 15600$ [\AA]', loc='center')

	# Use the hinode lines
	elif instrument == "Hinode":
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - 6300$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - 6300$ [\AA]', loc='center')
	else:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle1 = 0.51
		else:
			xtitle1 = 0.55
		if instrument == "GRIS":
			fig.suptitle("Near-Infrared Lines", y=0.98, x=xtitle1)
		elif instrument == "Hinode":
			fig.suptitle("Visible Lines", y=0.98, x=xtitle1)
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle1)
		elif title1 is not None:
			fig.suptitle(title1, y=0.98, x=xtitle1)

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
		plt.tight_layout(h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout()

	plt.savefig(savepath + "sir_inversion_stokes" + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-Bz"]
	index  = [1,2,3,4,5,6,7,8,9,10,11]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$","B [G]"]
	titles   = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength $B$", r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
			r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$", r"Magnetic Field Strength $B_z$"]

	i = 0

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8,7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)

			# Plot
			# Synthesised model = Real model
			if inputs[i] == "-Bz":
				ax1.plot(pars[0], pars[4]*np.cos(pars[6]/180*np.pi), label="Best Fit Model", color='#FF2C00')
				# Error of fit
				ax1.fill_between(pars[0], (pars[4] - pars_err[4])*np.cos(pars[6]/180*np.pi),
						(pars[4] + pars_err[4])*np.cos(pars[6]/180*np.pi), alpha = 0.5,
						color='#FF2C00', lw=0)
			else:
				ax1.plot(pars[0], pars[index[i]], label="Best Fit Model", color='#FF2C00')

				# Error of fit
				ax1.fill_between(pars[0], pars[index[i]] - pars_err[index[i]],
							pars[index[i]] + pars_err[index[i]], alpha = 0.5,
							color='#FF2C00', lw=0)

			# Set xlimits
			ax1.set_xlim(pars[0][0], pars[0][-1])

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{c}$")
			ax1.set_ylabel(labels[index[i]])

			if index[i] == 2 or index[i] == 9:
				ax1.semilogy()

			# Legend
			#ax1.legend()
		
			ax1.set_title(titles[i])
			# set the spacing between subplots
			plt.tight_layout(pad=2)
			plt.savefig(savepath + "sir_inversion_" + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = log_tau[-1]
	if lim_max < -3:
		lim_max = -3
		# Cut data so that the plot limits are adjusted to the shorted range
		I = np.where(log_tau < -3)[0][0]
		log_tau = log_tau[0:I]
		pars  = pars[:,0:I]
		pars_err  = pars_err[:,0:I]
		
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	llabel = "Best Fit"
	ax1.plot(log_tau, pars[1], label=llabel, color='#FF2C00')
	ax2.plot(log_tau, pars[4], label=llabel, color='#FF9500')
	ax3.plot(log_tau, pars[5], label=llabel, color='#474747')
	ax4.plot(log_tau, pars[6], label=llabel, color='#C20078')

	ax1.fill_between(pars[0], pars[1] - pars_err[1],
				pars[1] + pars_err[1], alpha = 0.5,
				color='#FF2C00', lw=0)
	ax2.fill_between(pars[0], pars[4] - pars_err[4],
				pars[4] + pars_err[4], alpha = 0.5,
				color='#FF9500', lw=0)
	ax3.fill_between(pars[0], pars[5] - pars_err[5],
				pars[5] + pars_err[5], alpha = 0.5,
				color='#474747', lw=0)
	ax4.fill_between(pars[0], pars[6] - pars_err[6],
				pars[6] + pars_err[6], alpha = 0.5,
				color='#C20078', lw=0)
	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(log_tau[0], lim_max)
	ax2.set_xlim(log_tau[0], lim_max)
	ax3.set_xlim(log_tau[0], lim_max)
	ax4.set_xlim(log_tau[0], lim_max)

	Min, Max = ax2.get_ylim()
	ax2.set_ylim(Min,Max*1.15)
	Min, Max = ax3.get_ylim()
	if abs(Min / Max) > 10:
		Max = Max / 1.15 * 3.2
	if Max < 0:
		Max = Max / 1.15 * 0.85
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

	#ax1.legend(loc='upper right')
	#ax2.legend(loc='upper right')
	#ax3.legend(loc='upper right')
	#ax4.legend(loc='upper right')

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

	
	plt.savefig(savepath + "sir_inversion_result" + add)


def inversion_single_2C(fit : str, obs : str, phy1 : str, phy2 : str):
	r"""
	
	Parameters
	----------
	fit : str
		Best Fit Profile
	obs : str
		Observation profiles
	phy1 : str
		Best Fit Model 1
	phy2 : str
		Best Fit Model 2

	Returns
	-------
	None

	Other Parameters
	----------------
	Additional parameters given as an argument when the script is executed.

	-save [str]
		Additional save path (optional, default './')
	-add [str]
			Additional text in filenames (optional)
	-label [str]
		Add label text (optional)
	-title [str]
		Title in Stokes plot
	-T
		Plot temperature in K
	-Pe
		Plot electron pressure in dyn/cm^2
	-vmicro
		Plot microturbulence in cm/s
	-B
		Plot magnetic field strength in Gauss
	-vlos
		Plot line of sight velocity in cm/s
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
	-syn
		Synthesised model .mod file
	-vertical
		Plot spectra vertically
	-num [int]
		Number of the line considered
	-line [str]
		Line file if -num is used
	-instrument [str]
		which instrument is used (GRIS 15648A line or Hinode) for relative wavelength, otherwise asked
	
	"""

	num = 0
	if '-num' in sys.argv:
		num = int(sys.argv[sys.argv.index("-num")+1])
	if '-line' in sys.argv:
		line = sys.argv[sys.argv.index("-line")+1]
	instrument =''
	if '-instrument' in sys.argv:
		instrument = sys.argv[sys.argv.index("-instrument")+1]

	# Observation from synthesis	
	ll, I, Q, U, V = sir.read_profile(obs, num)

	# Load fit data
	ll_fit, I_fit, Q_fit, U_fit, V_fit = sir.read_profile(fit, num)

	# Load physical parameters
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(phy1)
	pars1 = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(phy2)
	pars2 = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])

	# Load error of best fit model
	err1 = phy1.replace(".per",".err")
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(err1)
	pars1_err = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])
	err2 = phy2.replace(".per",".err")
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(err2)
	pars2_err = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])

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

	# Additional label
	add_label = '_'
	if '-label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]

	#############################################################
	#					PLOTTING STUFF					#
	#############################################################
	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Get colors used in the actual cycle

	# Change rel. to 6300 for Hinode lines
	if '-num' in sys.argv:
		Line = sir.read_line(line)
		Lines = Line['Line']
		ll0 = Line['wavelength']
		ll0 = ll0[np.where(Lines == num)[0][0]]

	elif instrument == "Hinode":
		print("The code assumes the line corresponds to 6301.5012 A! If you want to change it, take a look at line 124 in the script.")
		ll -= -6301.5012+6300.0000
		ll_fit -= -6301.5012+6300.0000
	elif instrument == "GRIS":
		print("The code assumes the line corresponds to 15648.515 A! If you want to change it, take a look at line 129 in the script.")
		ll -= -15648.514+15600.0000
		ll_fit -= -15648.514+15600.0000

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
	llabel = "Observation"
	ax1.plot(ll, I, "x", label=llabel)
	ax2.plot(ll, Q, "x", label=llabel)
	ax3.plot(ll, U, "x", label=llabel)
	ax4.plot(ll, V, "x", label=llabel)
	ax1.plot(ll, I, "-", color=colors[0], alpha = 0.5)
	ax2.plot(ll, Q, "-", color=colors[0], alpha = 0.5)
	ax3.plot(ll, U, "-", color=colors[0], alpha = 0.5)
	ax4.plot(ll, V, "-", color=colors[0], alpha = 0.5)

	ax1.plot(ll, I_fit, "-", label = "Best Fit", color=colors[1])
	ax2.plot(ll, Q_fit, "-", label = "Best Fit", color=colors[1])
	ax3.plot(ll, U_fit, "-", label = "Best Fit", color=colors[1])
	ax4.plot(ll, V_fit, "-", label = "Best Fit", color=colors[1])

	# Set xlimits
	ax1.set_xlim(ll[0], ll[-1])
	ax2.set_xlim(ll[0], ll[-1])
	ax3.set_xlim(ll[0], ll[-1])
	ax4.set_xlim(ll[0], ll[-1])

	#######################################################################
	# Set labels												#
	# The labels depend on the arguments and the chosen line			#
	# The code is so long as I tried to make it as automatic as possible	#
	#######################################################################
	ax1.set_ylabel(r'$I / I_c$')
	ax2.set_ylabel(r'$Q / I_c$')
	ax3.set_ylabel(r'$U / I_c$')
	ax4.set_ylabel(r'$V / I_c$')
	# If -line and -num in arguments => determine rel. wavelength from there
	# Do not include if hinode lines are selected because then the rel. is 6300
	# to make it easier to be read from the plot
	title1 = None
	if '-num' in sys.argv:
		Line = sir.read_line(line)
		Lines = Line['Line']
		ll0 = Line['wavelength']
		ll0 = ll0[np.where(Lines == num)[0][0]]
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + str(ll0) + r' [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - $' + str(ll0) + r' [\AA]', loc='center')
		title1 = "Fe I " + str(ll0)

	# Use the four 1.5 lines
	elif instrument == "GRIS":
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - 15600$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - 15600$ [\AA]', loc='center')

	# Use the hinode lines
	elif instrument == "Hinode":
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - 6300$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - 6300$ [\AA]', loc='center')
	else:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle1 = 0.51
		else:
			xtitle1 = 0.55
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle1)
		elif title1 is not None:
			fig.suptitle(title1, y=0.98, x=xtitle1)

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
		plt.tight_layout(h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout()

	plt.savefig(savepath + "inversion_stokes" + add)

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

			# Plot
			# Synthesised model = Real model
			ax1.plot(pars1[0], pars1[index[i]], label="Best Fit Model 1", color=colors[0])
			ax1.plot(pars2[0], pars2[index[i]], label="Best Fit Model 2", color=colors[1])

			# Error of fit
			ax1.fill_between(pars1[0], pars1[index[i]] - pars1_err[index[i]],
						pars1[index[i]] + pars1_err[index[i]], alpha = 0.5,
						color=colors[0], lw=0)
			ax1.fill_between(pars2[0], pars2[index[i]] - pars2_err[index[i]],
						pars2[index[i]] + pars2_err[index[i]], alpha = 0.5,
						color=colors[1], lw=0)

			# Set xlimits
			ax1.set_xlim(pars1[0][0], pars1[0][-1])

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
			plt.savefig(savepath + "inversion_" + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = pars1[0][-1]
	if lim_max < -3:
		lim_max = -3
		# Cut data so that the plot limits are adjusted to the shorted range
		I = np.where(pars1[0] < -3)[0][0]
		print(pars1.shape)
		#pars1[0] = pars1[0][0:I]
		pars1  = pars1[:,0:I]
		pars1_err  = pars1_err[:,0:I]
		I = np.where(pars2[0] < -3)[0][0]
		#pars2[0] = pars2[0][0:I]
		pars2  = pars2[:,0:I]
		pars2_err  = pars2_err[:,0:I]


	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	llabel = "Best Fit Model 1"
	ax1.plot(pars1[0], pars1[1], label=llabel, color=colors[0])
	ax2.plot(pars1[0], pars1[4], label=llabel, color=colors[0])
	ax3.plot(pars1[0], pars1[5], label=llabel, color=colors[0])
	ax4.plot(pars1[0], pars1[6], label=llabel, color=colors[0])
	llabel = "Best Fit Model 2"
	ax1.plot(pars2[0], pars2[1], label=llabel, color=colors[1])
	ax2.plot(pars2[0], pars2[4], label=llabel, color=colors[1])
	ax3.plot(pars2[0], pars2[5], label=llabel, color=colors[1])
	ax4.plot(pars2[0], pars2[6], label=llabel, color=colors[1])

	ax1.fill_between(pars1[0], pars1[1] - pars1_err[1],
				pars1[1] + pars1_err[1], alpha = 0.5,
				color=colors[0], lw=0)
	ax2.fill_between(pars1[0], pars1[4] - pars1_err[4],
				pars1[4] + pars1_err[4], alpha = 0.5,
				color=colors[0], lw=0)
	ax3.fill_between(pars1[0], pars1[5] - pars1_err[5],
				pars1[5] + pars1_err[5], alpha = 0.5,
				color=colors[0], lw=0)
	ax4.fill_between(pars1[0], pars1[6] - pars1_err[6],
				pars1[6] + pars1_err[6], alpha = 0.5,
				color=colors[0], lw=0)
	ax1.fill_between(pars2[0], pars2[1] - pars2_err[1],
				pars2[1] + pars2_err[1], alpha = 0.5,
				color=colors[1], lw=0)
	ax2.fill_between(pars2[0], pars2[4] - pars2_err[4],
				pars2[4] + pars2_err[4], alpha = 0.5,
				color=colors[1], lw=0)
	ax3.fill_between(pars2[0], pars2[5] - pars2_err[5],
				pars2[5] + pars2_err[5], alpha = 0.5,
				color=colors[1], lw=0)
	ax4.fill_between(pars2[0], pars2[6] - pars2_err[6],
				pars2[6] + pars2_err[6], alpha = 0.5,
				color=colors[1], lw=0)
	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(pars1[0][0], lim_max)
	ax2.set_xlim(pars1[0][0], lim_max)
	ax3.set_xlim(pars1[0][0], lim_max)
	ax4.set_xlim(pars1[0][0], lim_max)

	Min, Max = ax2.get_ylim()
	ax2.set_ylim(Min,Max*1.15)
	Min, Max = ax3.get_ylim()
	if abs(Min / Max) > 10:
		Max = Max / 1.15 * 3.2
	if Max < 0:
		Max = Max / 1.15 * 0.85
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
	#ax2.legend(loc='upper right')
	#ax3.legend(loc='upper right')
	#ax4.legend(loc='upper right')

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


	plt.savefig(savepath + "inversion_result" + add)






def inversion_single_mc(fit : str, obs : str, phy : str, syn : str):
	"""
	Print a single inversion directly from the .mod and .per files

	Parameters
	----------
	fit : str
		Best Fit Profile
	obs : str
		Observation profiles
	phy : str
		Best Fit Model
	syn : str
		Synthesis Model

	Returns
	-------
	None

	"""
	# Import library
	sir.mpl_library()

	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	

	num = 0
	if '-num' in sys.argv:
		num = int(sys.argv[sys.argv.index("-num")+1])
	if '-line' in sys.argv:
		line = sys.argv[sys.argv.index("-line")+1]


	# Observation from synthesis	
	ll, I, Q, U, V = sir.read_profile(obs, num)

	# Load fit data
	ll_fit, I_fit, Q_fit, U_fit, V_fit = sir.read_profile(fit, num)

	# Load physical parameters
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(phy)
	pars = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])

	# Load physical parameter of synthesised model
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(syn)
	syn = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])

	# Load error of best fit model
	err = fit.replace(".per",".err")
	log_tau, T, Pe, v_micro, B, vlos, inc, azimuth, z, Pg, rho = sir.read_model(err)
	pars_err = np.array([log_tau, T, Pe, v_micro, B, vlos/1e5, inc, azimuth, z, Pg, rho])


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

	# Additional label
	add_label = '_'
	if '-label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]

	#####################################################
	#					PLOTTING STUFF					#
	#####################################################

	# Change rel. to 6300 for Hinode lines
	if '-num' in sys.argv:
		Line = sir.read_line(line)
		Lines = Line['Line']
		ll0 = Line['wavelength']
		ll0 = ll0[np.where(Lines == num)[0][0]]

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
	llabel = "Synthesised Profile"
	ax1.plot(ll, I, "x", label=llabel)
	ax2.plot(ll, Q, "x", label=llabel)
	ax3.plot(ll, U, "x", label=llabel)
	ax4.plot(ll, V, "x", label=llabel)

	ax1.plot(ll, I_fit, "-", label = "Best Fit")
	ax2.plot(ll, Q_fit, "-", label = "Best Fit")
	ax3.plot(ll, U_fit, "-", label = "Best Fit")
	ax4.plot(ll, V_fit, "-", label = "Best Fit")

	# Set xlimits
	ax1.set_xlim(ll[0], ll[-1])
	ax2.set_xlim(ll[0], ll[-1])
	ax3.set_xlim(ll[0], ll[-1])
	ax4.set_xlim(ll[0], ll[-1])

	#######################################################################
	# Set labels												#
	# The labels depend on the arguments and the chosen line			#
	# The code is so long as I tried to make it as automatic as possible	#
	#######################################################################
	ax1.set_ylabel(r'$I / I_c$')
	ax2.set_ylabel(r'$Q / I_c$')
	ax3.set_ylabel(r'$U / I_c$')
	ax4.set_ylabel(r'$V / I_c$')
	# If -line and -num in arguments => determine rel. wavelength from there
	# Do not include if hinode lines are selected because then the rel. is 6300
	# to make it easier to be read from the plot
	title1 = None
	if '-num' in sys.argv:
		Line = sir.read_line(line)
		Lines = Line['Line']
		ll0 = Line['wavelength']
		ll0 = ll0[np.where(Lines == num)[0][0]]
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + str(ll0) + r' [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - $' + str(ll0) + r' [\AA]', loc='center')
		title1 = "Fe I " + str(ll0)

	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')
	ax4.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle1 = 0.51
		else:
			xtitle1 = 0.55
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle1)
		elif title1 is not None:
			fig.suptitle(title1, y=0.98, x=xtitle1)

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
		plt.tight_layout(h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout()

	plt.savefig(savepath + "sir_inversion_stokes" + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-Bz"]
	index  = [1,2,3,4,5,6,7,8,9,10,11]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$","B [G]"]
	titles   = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength $B$", r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
			r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$", r"Magnetic Field Strength $B_z$"]

	i = 0

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8,7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)

			# Plot
			# Synthesised model = Real model
			if inputs[i] == "-Bz":
				ax1.plot(syn[0], syn[4]*np.cos(syn[6]/180*np.pi), label="Model")
				ax1.plot(pars[0], pars[4]*np.cos(pars[6]/180*np.pi), label="Best Fit Model", color='#FF2C00')
				# Error of fit
				ax1.fill_between(pars[0], (pars[4] - pars_err[4])*np.cos(pars[6]/180*np.pi),
						(pars[4] + pars_err[4])*np.cos(pars[6]/180*np.pi), alpha = 0.5,
						color='#FF2C00', lw=0)
			else:
				ax1.plot(syn[0], syn[index[i]], label="Model")
				ax1.plot(pars[0], pars[index[i]], label="Best Fit Model", color='#FF2C00')

				# Error of fit
				ax1.fill_between(pars[0], pars[index[i]] - pars_err[index[i]],
							pars[index[i]] + pars_err[index[i]], alpha = 0.5,
							color='#FF2C00', lw=0)

			# Set xlimits
			ax1.set_xlim(pars[0][0], pars[0][-1])

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
			plt.savefig(savepath + "sir_inversion_" + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = log_tau[-1]
	if lim_max < -3:
		lim_max = -3
		# Cut data so that the plot limits are adjusted to the shorted range
		I = np.where(log_tau < -3)[0][0]
		log_tau = log_tau[0:I]
		pars  = pars[:,0:I]
		pars_err  = pars_err[:,0:I]
		syn  = syn [:,0:I]
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	llabel = "Synthesised Model"
	ax1.plot(log_tau, syn[1], label=llabel, color='#0C5DA5')
	ax2.plot(log_tau, syn[4], label=llabel, color='#00B945')
	ax3.plot(log_tau, syn[5], label=llabel, color='#845B97')
	ax4.plot(log_tau, syn[6], label=llabel, color='#054907')

	llabel = "Best Fit"
	ax1.plot(log_tau, pars[1], label=llabel, color='#FF2C00')
	ax2.plot(log_tau, pars[4], label=llabel, color='#FF9500')
	ax3.plot(log_tau, pars[5], label=llabel, color='#474747')
	ax4.plot(log_tau, pars[6], label=llabel, color='#C20078')

	ax1.fill_between(pars[0], pars[1] - pars_err[1],
				pars[1] + pars_err[1], alpha = 0.5,
				color='#FF2C00', lw=0)
	ax2.fill_between(pars[0], pars[4] - pars_err[4],
				pars[4] + pars_err[4], alpha = 0.5,
				color='#FF9500', lw=0)
	ax3.fill_between(pars[0], pars[5] - pars_err[5],
				pars[5] + pars_err[5], alpha = 0.5,
				color='#474747', lw=0)
	ax4.fill_between(pars[0], pars[6] - pars_err[6],
				pars[6] + pars_err[6], alpha = 0.5,
				color='#C20078', lw=0)
	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(log_tau[0], lim_max)
	ax2.set_xlim(log_tau[0], lim_max)
	ax3.set_xlim(log_tau[0], lim_max)
	ax4.set_xlim(log_tau[0], lim_max)

	Min, Max = ax2.get_ylim()
	ax2.set_ylim(Min,Max*1.15)
	Min, Max = ax3.get_ylim()
	if abs(Min / Max) > 10:
		Max = Max / 1.15 * 3.2
	if Max < 0:
		Max = Max / 1.15 * 0.85
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
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle2)

	if "-vertical" in sys.argv:	
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		# set the spacing between subplots
		plt.tight_layout(pad=2)

	
	plt.savefig(savepath + "sir_inversion_result" + add)



# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()

	if sys.argv[1] == '1C':
		fit = sys.argv[2]   # Best fit
		obs = sys.argv[3]   # Synthesis profile
		phy = sys.argv[4]   # Physical parameter of best fit
		inversion_single_1C(fit, obs, phy)

	elif sys.argv[1] == '2C':
		fit  = sys.argv[2]   # Best fit
		obs  = sys.argv[3]   # Observed profile
		phy1 = sys.argv[4]   # Physical parameter of best fit
		phy2 = sys.argv[5]   # Physical parameter of best fit
		inversion_single_2C(fit, obs, phy1, phy2)

	elif sys.argv[1] == "MC":
		fit = sys.argv[2]   # Best fit
		obs = sys.argv[3]   # Synthesis profile
		phy = sys.argv[4]   # Physical parameter of best fit
		syn = sys.argv[5]   # Synthesis model

		inversion_single_mc(fit, obs, phy, syn)
	else:
		print(f"[result] Mode '{sys.argv[1]}' not known or undefined!")

	