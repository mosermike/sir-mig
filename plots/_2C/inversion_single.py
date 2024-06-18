"""
Plots the result of the SIR synthesis
"""

import numpy as np 
import sys, os
sys.path.append(sys.path[0] + "/../../src")
import sir
import definitions as d
import matplotlib.pyplot as plt
from os.path import exists

# Import library
dirname = os.path.split(os.path.abspath(__file__))[0]
plt.rcParams["savefig.format"] = "pdf"
if d.plt_lib != "":
	plt.style.use(d.plt_lib)
else:
	if exists(dirname + '/../mml.mplstyle'):
		plt.style.use(dirname + '/../mml.mplstyle')
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


if "-h" in sys.argv or (len(sys.argv) < 5):
	print("plot_inversion_single - Plots the result of a inversion without using a config file")
	print("Usage: python plot_inversion_single [OPTION]")
	print()
	print("1:         Observation Profile")
	print("2:         Best Fit Profile")
	print("3:         Best Fit Model 1")
	print("4:         Best Fit Model 2")

	print("-save:     Additional save path (optional, default './')")
	print("-add:      Additional text in filenames (optional)")
	print("-label:    Add label text (optional)")
	print("-title:    Title in Stokes plot")
	print("-T:        Plot temperature in K")
	print("-Pe:       Plot electron pressure in dyn/cm^2")
	print("-vmicro:   Plot microturbulence in cm/s")
	print("-B:        Plot magentic field strength in Gauss")
	print("-vlos:     Plot line of sight velocity in cm/s")
	print("-inc:      Plot inclination by subtracting in deg")
	print("-azi:      Plot azimuth by adding in deg")
	print("-z:        Plot real height in km")
	print("-Pg:       Plot gas pressure in dyn/cm^2")
	print("-rho:      Plot density")
	print("-syn:      Synthesised model .mod file")
	print("-vertical: Plot spectra vertically")
	print("-num:      Number of the line considered")
	print("-line:     Line file if -num is used")
	if "-h" not in sys.argv:
		print("[ERROR] Not enough arguments passed!")
	sys.exit()

#############################################################
#			READ INPUT AND LOAD DATA					#
#############################################################
obs  = sys.argv[2]   # Observed profile
fit  = sys.argv[1]   # Best fit
phy1 = sys.argv[3]   # Physical parameter of best fit
phy2 = sys.argv[4]   # Physical parameter of best fit

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

  
plt.savefig(savepath + "inversion_result" + add)





