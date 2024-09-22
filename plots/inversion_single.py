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
	sir.option("[5. Pos]","Synthesis Model or Best Fit Model 2C (Modes: MC, 2C)")
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
	sir.option("-line","Line file for using absolute wavelength values")
	sir.option("-instrument [str]","which instrument is used for relative wavelength (Hinode or GRIS at 15648A) (Mode `1C` and `2C`), otherwise asked")
	sir.option("-err","Print errorbars")
	sys.exit()


def inversion_single_1C(fit_str : str, obs_str : str, phy_str : str, syn_str : str=""):
	r"""
	
	Parameters
	----------
	fit_str : str
		Best Fit Profile
	obs_str : str
		Observation profiles
	phy_str : str
		Best Fit Model
	syn_str : str, optional
		Used for the synthesis model, e.g. in mode MC, by default "".

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
		Line file for using absolute wavelength values
	
	
	"""

	# Import library
	sir.mpl_library()
	

	import profile_stk as p
	# Observation from synthesis	
	obs = p.read_profile(obs_str)

	# Load fit data
	fit = p.read_profile(fit_str)

	import model_atm as m
	# Load physical parameters
	phy = m.read_model(phy_str)

	if syn_str != "":
		syn = m.read_model(syn_str)

	# Load error of best fit model
	err = m.read_model(phy_str.replace(".per",".err"))
	
	if '-line' in sys.argv:
		line = sys.argv[sys.argv.index("-line")+1]
		# Determine line numbers
		nums = []
		for i in obs.indx:
			if i not in nums:
				nums.append(i)
		
		# Add the line core to the wavelenght numbers
		for i in nums:
			ll0 = sir.determine_line_core(line, i)
			obs.wave[obs.indx == i] += ll0
			fit.wave[fit.indx == i] += ll0
	
	
	# Change to abs. wavelength to the line core of the first number
	if "-num" in sys.argv:
		num = int(sys.argv[sys.argv.index("-num")+1])
		obs.cut_to_wave(np.array([[obs.wave[obs.indx == num][0], obs.wave[obs.indx == num][1]-obs.wave[obs.indx == num][0],len(obs.wave[obs.indx == num])]]))
		fit.cut_to_wave(np.array([[fit.wave[fit.indx == num][0], fit.wave[fit.indx == num][1]-fit.wave[fit.indx == num][0],len(fit.wave[fit.indx == num])]]))
	
	label_x = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if label_x != '0':
		obs.wave -= float(label_x)
		fit.wave -= float(label_x)

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

	########################
	#  Plot I, Q, U and V  #
	########################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,14), sharex=True,
			gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0, layout="compressed")
	elif "-hor" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, figsize=(17.39,4.31), layout="compressed")
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), layout="compressed")


	# Add label
	ax2.scatter([], [], color="w", alpha=0, label=add_label)

	############################
	# Plot the Stokes profiles #
	############################
	llabel = "Observation"
	ax1.plot(obs.wave, obs.stki[0,0], "x", label=llabel)
	ax2.plot(obs.wave, obs.stkq[0,0], "x", label=llabel)
	ax3.plot(obs.wave, obs.stku[0,0], "x", label=llabel)
	ax4.plot(obs.wave, obs.stkv[0,0], "x", label=llabel)

	ax1.plot(fit.wave, fit.stki[0,0], "-", label = "Best Fit")
	ax2.plot(fit.wave, fit.stkq[0,0], "-", label = "Best Fit")
	ax3.plot(fit.wave, fit.stku[0,0], "-", label = "Best Fit")
	ax4.plot(fit.wave, fit.stkv[0,0], "-", label = "Best Fit")

	# Set xlimits
	ax1.set_xlim(fit.wave[0], fit.wave[-1])
	ax2.set_xlim(fit.wave[0], fit.wave[-1])
	ax3.set_xlim(fit.wave[0], fit.wave[-1])
	ax4.set_xlim(fit.wave[0], fit.wave[-1])

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
	if label_x != "0":
		if "-hor" in sys.argv:
			ax1.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA') 
			ax2.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
		ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
	else:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\lambda$' + r' [\AA]')
		ax4.set_xlabel(r'$\lambda$' + r' [\AA]')

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
		ax1.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots
	else:
		ax1.legend()
		# set the spacing between subplots	
		

	fig.savefig(savepath + "inversion_stokes" + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-Bz"]
	index  = [1,2,3,4,5,6,7,8,9,10,11]
	labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$","B [G]"]
	titles   = [r"Temperature T", r"Electron Pressure $P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength $B$", r"LOS Velocity $\mathrm{v}_{\mathrm{los}}$",
			r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $P_g$", r"Density $\rho$", r"Magnetic Field Strength $B_z$"]

	i = 0
	pars = [phy.tau,phy.T,phy.Pe,phy.vmicro,phy.B,phy.vlos,phy.gamma,phy.phi,phy.z,phy.Pg,phy.rho]
	pars_err = [err.tau,err.T,err.Pe,err.vmicro,err.B,err.vlos,err.gamma,err.phi,err.z,err.Pg,err.rho]

	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8,7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)

			# Plot
			# Synthesised model = Real model
			if inputs[i] == "-Bz":
				ax1.plot(phy.tau, phy.B[0,0]*np.cos(phy.gamma[0,0]/180*np.pi), label="Best Fit Model", color='#FF2C00')
				# Error of fit
				if "-err" in sys.argv:
					ax1.fill_between(pars[0], (phy.B[0,0] - err.B[0,0])*np.cos(phy.gamma[0,0]/180*np.pi),
						(phy.B[0,0] + err.B[0,0])*np.cos(phy.gamma[0,0]/180*np.pi), alpha = 0.5,
						color='#FF2C00', lw=0)
			else:
				ax1.plot(pars[0], pars[index[i]][0,0], label="Best Fit Model", color='#FF2C00')

				# Error of fit
				if "-err" in sys.argv:
					ax1.fill_between(pars[0], pars[index[i]][0,0] - pars_err[index[i]][0,0],
							pars[index[i]][0,0] + pars_err[index[i]][0,0], alpha = 0.5,
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
			plt.savefig(savepath + "inversion_" + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	
	lim_max = -3.0
	phy.set_limit(lim_max)
	err.set_limit(lim_max)
	
	pars = [phy.tau,phy.T,phy.Pe,phy.vmicro,phy.B,phy.vlos,phy.gamma,phy.phi,phy.z,phy.Pg,phy.rho]
	pars_err = [err.tau,err.T,err.Pe,err.vmicro,err.B,err.vlos,err.gamma,err.phi,err.z,err.Pg,err.rho]
	
	colors = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Get colors used in the actual cycle

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			gridspec_kw=dict(hspace=0), layout="compressed")
		fig.subplots_adjust(hspace=0, wspace=0)
	elif "-hor" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, figsize=(17.39,4.31), layout="compressed")
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True, layout="compressed")

	inds = [0,1,2,3]
	if syn_str != "":
		inds = [1,1,1,1]
		llabel = "Synthesis"
		ax1.plot(syn.tau, syn.T[0,0], label=llabel, color=colors[0])
		ax2.plot(syn.tau, syn.B[0,0], label=llabel, color=colors[0])
		ax3.plot(syn.tau, syn.vlos[0,0], label=llabel, color=colors[0])
		ax4.plot(syn.tau, syn.gamma[0,0], label=llabel, color=colors[0])

	llabel = "Best Fit"
	ax1.plot(phy.tau, phy.T[0,0], label=llabel, color=colors[inds[0]])
	ax2.plot(phy.tau, phy.B[0,0], label=llabel, color=colors[inds[1]])
	ax3.plot(phy.tau, phy.vlos[0,0], label=llabel, color=colors[inds[2]])
	ax4.plot(phy.tau, phy.gamma[0,0], label=llabel, color=colors[inds[3]])

	
	if "-err" in sys.argv:
		ax1.fill_between(pars[0], pars[1][0,0] - pars_err[1][0,0],
					pars[1][0,0] + pars_err[1][0,0], alpha = 0.5,
					color=colors[0], lw=0)
		ax2.fill_between(pars[0], pars[4][0,0] - pars_err[4][0,0],
					pars[4][0,0] + pars_err[4][0,0], alpha = 0.5,
					color=colors[1], lw=0)
		ax3.fill_between(pars[0], pars[5][0,0] - pars_err[5][0,0],
					pars[5][0,0] + pars_err[5][0,0], alpha = 0.5,
					color=colors[2], lw=0)
		ax4.fill_between(pars[0], pars[6][0,0] - pars_err[6][0,0],
					pars[6][0,0] + pars_err[6][0,0], alpha = 0.5,
					color=colors[3], lw=0)
		#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(phy.tau[0], lim_max)
	ax2.set_xlim(phy.tau[0], lim_max)
	ax3.set_xlim(phy.tau[0], lim_max)
	ax4.set_xlim(phy.tau[0], lim_max)

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
	if "-hor" not in sys.argv:
		ax1.set_xlabel(r"$\log \tau_{c}$")
		ax2.set_xlabel(r"$\log \tau_{c}$")
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

	if syn_str != "":
		ax1.legend()
	

	# Set title position depending on the chosen plot and consider the flags hinode and gris
	if title != "-1":
		if "-vertical" in sys.argv:
			xtitle2 = 0.51
		else:
			xtitle2 = 0.55
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle2)


	
	fig.savefig(savepath + "inversion_result" + add)


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
		Line file for using absolute wavelength values
	
	
	"""
	sir.mpl_library()
	

	import profile_stk as p
	# Observation from synthesis	
	obs1 = p.read_profile(obs)

	# Load fit data
	fit1 = p.read_profile(fit)

	import model_atm as m
	# Load physical parameters
	phy11 = m.read_model(phy1)
	phy22 = m.read_model(phy2)

	# Load error of best fit model
	err1 = m.read_model(phy1.replace(".per",".err"))
	err2 = m.read_model(phy2.replace(".per",".err"))
	
	if '-line' in sys.argv:
		line = sys.argv[sys.argv.index("-line")+1]
		# Determine line numbers
		nums = []
		for i in obs1.indx:
			if i not in nums:
				nums.append(i)
		
		# Add the line core to the wavelenght numbers
		for i in nums:
			ll0 = sir.determine_line_core(line, i)
			obs1.wave[obs1.indx == i] += ll0
			fit1.wave[fit1.indx == i] += ll0


	# Change to abs. wavelength to the line core of the first number
	if "-num" in sys.argv:
		
		num = int(sys.argv[sys.argv.index("-num")+1])
		obs1.cut_to_wave([obs1.wave[obs1.indx == num][0], obs1.wave[obs1.indx == num][1]-obs1.wave[obs1.indx == num][0],len(obs1.wave[obs1.indx == num])])
		fit1.cut_to_wave([fit1.wave[fit1.indx == num][0], fit1.wave[fit1.indx == num][1]-fit1.wave[fit1.indx == num][0],len(fit1.wave[fit1.indx == num])])

	label_x = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if label_x != '0':
		obs1.wave -= float(label_x)
		fit1.wave -= float(label_x)
	
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

	########################
	#  Plot I, Q, U and V  #
	########################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,14), sharex=True,
			gridspec_kw=dict(hspace=0), layout="compressed")
		fig.subplots_adjust(hspace=0, wspace=0)
	elif "-hor" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, figsize=(17.39,4.31), layout="compressed")
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), layout="compressed")


	# Add label
	ax2.scatter([], [], color="w", alpha=0, label=add_label)

	############################
	# Plot the Stokes profiles #
	############################
	llabel = "Observation"
	ax1.plot(obs1.wave, obs1.stki, "x", label=llabel)
	ax2.plot(obs1.wave, obs1.stkq, "x", label=llabel)
	ax3.plot(obs1.wave, obs1.stku, "x", label=llabel)
	ax4.plot(obs1.wave, obs1.stkv, "x", label=llabel)
	ax1.plot(obs1.wave, obs1.stki, "-", color=colors[0], alpha = 0.5)
	ax2.plot(obs1.wave, obs1.stkq, "-", color=colors[0], alpha = 0.5)
	ax3.plot(obs1.wave, obs1.stku, "-", color=colors[0], alpha = 0.5)
	ax4.plot(obs1.wave, obs1.stkv, "-", color=colors[0], alpha = 0.5)

	ax1.plot(fit1.wave, fit1.stki, "-", label = "Best Fit", color=colors[1])
	ax2.plot(fit1.wave, fit1.stkq, "-", label = "Best Fit", color=colors[1])
	ax3.plot(fit1.wave, fit1.stku, "-", label = "Best Fit", color=colors[1])
	ax4.plot(fit1.wave, fit1.stkv, "-", label = "Best Fit", color=colors[1])

	# Set xlimits
	ax1.set_xlim(obs1.wave[0], obs1.wave[-1])
	ax2.set_xlim(obs1.wave[0], obs1.wave[-1])
	ax3.set_xlim(obs1.wave[0], obs1.wave[-1])
	ax4.set_xlim(obs1.wave[0], obs1.wave[-1])

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
	if label_x != "0":
		if "-hor" in sys.argv:
			ax1.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
			ax2.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
		ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')
	else:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\lambda$' + r' [\AA]')
		ax4.set_xlabel(r'$\lambda$' + r' [\AA]')

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
		ax1.legend(bbox_to_anchor=(1.01,0.95))
	elif "-hor" in sys.argv:
		ax1.legend()
	else:
		ax1.legend()
		# set the spacing between subplots	
		

	fig.savefig(savepath + "inversion_stokes" + add)

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

	phy = phy11.copy()
	err = err1.copy()
	pars1 = [phy.tau[0,0],phy.T[0,0],phy.Pe[0,0],phy.vmic[0,0],phy.B[0,0],phy.vlos[0,0],phy.gamma[0,0],phy.phi[0,0],phy.z[0,0],phy.Pg[0,0],phy.rho[0,0]]
	pars1_err = [err.tau[0,0],err.T[0,0],err.Pe[0,0],err.vmic[0,0],err.B[0,0],err.vlos[0,0],err.gamma[0,0],err.phi[0,0],err.z[0,0],err.Pg[0,0],err.rho[0,0]]

	phy = phy22.copy()
	err = err2.copy()
	pars2 = [phy.tau[0,0],phy.T[0,0],phy.Pe[0,0],phy.vmic[0,0],phy.B[0,0],phy.vlos[0,0],phy.gamma[0,0],phy.phi[0,0],phy.z[0,0],phy.Pg[0,0],phy.rho[0,0]]
	pars2_err = [err.tau[0,0],err.T[0,0],err.Pe[0,0],err.vmic[0,0],err.B[0,0],err.vlos[0,0],err.gamma[0,0],err.phi[0,0],err.z[0,0],err.Pg[0,0],err.rho[0,0]]

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
			if "-err" in sys.argv:
				ax1.fill_between(pars1[0], pars1[index[i]] - pars1_err[index[i]],
							pars1[index[i]] + pars1_err[index[i]], alpha = 0.5,
							color=colors[0], lw=0)
				ax1.fill_between(pars2[0], pars2[index[i]] - pars2_err[index[i]],
							pars2[index[i]] + pars2_err[index[i]], alpha = 0.5,
							color=colors[1], lw=0)

			# Set xlimits
			ax1.set_xlim(pars1[0][0], pars1[0][-1])

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
			plt.savefig(savepath + "inversion_" + str(inputs[i][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max  = -3
	phy = phy11.copy()
	err = err1.copy()
	phy.set_lim(lim_max)
	err.set_lim(lim_max)
	pars1 = [phy.tau[0,0],phy.T[0,0],phy.Pe[0,0],phy.vmic[0,0],phy.B[0,0],phy.vlos[0,0],phy.gamma[0,0],phy.phi[0,0],phy.z[0,0],phy.Pg[0,0],phy.rho[0,0]]
	pars1_err = [err.tau[0,0],err.T[0,0],err.Pe[0,0],err.vmic[0,0],err.B[0,0],err.vlos[0,0],err.gamma[0,0],err.phi[0,0],err.z[0,0],err.Pg[0,0],err.rho[0,0]]

	phy = phy22.copy()
	err = err2.copy()
	phy.set_lim(lim_max)
	err.set_lim(lim_max)
	pars2 = [phy.tau[0,0],phy.T[0,0],phy.Pe[0,0],phy.vmic[0,0],phy.B[0,0],phy.vlos[0,0],phy.gamma[0,0],phy.phi[0,0],phy.z[0,0],phy.Pg[0,0],phy.rho[0,0]]
	pars2_err = [err.tau[0,0],err.T[0,0],err.Pe[0,0],err.vmic[0,0],err.B[0,0],err.vlos[0,0],err.gamma[0,0],err.phi[0,0],err.z[0,0],err.Pg[0,0],err.rho[0,0]]


	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0, layout="compressed")
	elif "-hor" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, figsize=(17.39,4.31), layout="compressed")
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True, layout="compressed")

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

	if "-err" in sys.argv:
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
	if "-hor" in sys.argv:
		ax1.set_xlabel(r'$\log\tau_c$')
		ax2.set_xlabel(r'$\log\tau_c$')
	if "-vertical" not in sys.argv:
		ax3.set_xlabel(r"$\log \tau_c$")
	ax4.set_xlabel(r"$\log \tau_c$")

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


	fig.savefig(savepath + "inversion_result" + add)





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

		inversion_single_1C(fit, obs, phy, syn)
	else:
		print(f"[result] Mode '{sys.argv[1]}' not known or undefined!")

	