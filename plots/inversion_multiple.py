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
	print("plot_inversion_multiple - Plots the result of multiple inversions")
	print("Usage: python plot_inversion_2 [OPTION]")
	print()
	sir.option("[1. Pos]","Config Models as a string config1.txt,config2.txt,...")
	sir.option("[2. Pos]","x positions as a string x1,x2,x3,x4,...")
	sir.option("[3. Pos]","y positions as a string y1,y2,y3,... (choose 0 for mode 'MC' and 'SY')")
	
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-add_label","Additional label in plots")
	sir.option("-labels", "Labels as a string where the string is split at ','.")
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

def inversion_multiple(confs : list, x : list, y : list, labels : list):
	""""
	Plots the result of multiple inversions. The code assumes the same x and y position for all inversions. In addition, it also assumes that in mode "MC"
	the synthesis model is the same

	Parameters
	----------
	confs : list
		List with dictionaries of the config files
	x : list
		x positions
	y : list
		y positions
	labels : list
		List of strings for the labels
	
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
	-vertical
		Plot it vertically
	-num [int]
		which line number is used (only for mode `MC`)
	-one
		Observations are the same and therefore only plotted once
	-err
		Plot errorbars

	Raises
	------
	RuntimeError
		if the lenghts of confs, x and y are not the same
	ValueError
		if the mode in the configs are not the same
	"""
	# Import library
	sir.mpl_library()

	if(len(x) != len(y) != len(confs)):
		raise RuntimeError("[inversion_multiple] Lengths of the arrays confs, x and y not the same!")
	
	if(np.all(np.array([i["mode"] for i in confs], dtype=str) != confs[0]["mode"])):
		raise ValueError(f"The mode must be the same for all configs. The first one is '{confs[0]['mode']}'.")
	
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################

	if confs[0]['mode'] == "MC":
		syns = [m.read_model(os.path.join(confs[i]["path"],confs[i]["syn_out"] + d.end_models)) for i in range(len(confs))]	# Synthesis Models 1
	elif confs[0]["mode"] == "SY":
		syns = [m.read_model(os.path.join(confs[i]["path"],confs[i]["syn_in"])) for i in range(len(confs))]
	
	if confs[0]['mode'] == "2C":
		syns = [m.read_model(os.path.join(confs[i]["path"],confs[i]["inv_out"] + d.end_models1)) for i in range(len(confs))]	# Synthesis Models 2	
		phys = [m.read_model(os.path.join(confs[i]["path"],confs[i]["inv_out"] + d.end_models2)) for i in range(len(confs))]	# Inversion Models 2
	elif confs[0]["mode"] != "SY":
		phys = [m.read_model(os.path.join(confs[i]["path"],confs[i]["inv_out"] + d.end_models)) for i in range(len(confs))]	# Inversion Models 1
		
	if confs[0]['mode'] == "MC":
		obss = [p.read_profile(os.path.join(confs[i]["path"],confs[i]["syn_out"] + d.end_stokes)) for i in range(len(confs))]	# Synthesis Profiles 1
	elif confs[0]["mode"] != "SY":
		obss = [p.read_profile(os.path.join(confs[i]["path"],confs["cube"])) for i in range(len(confs))]

	if confs[0]["mode"] == "SY":
		fits = [p.read_profile(os.path.join(confs[i]["path"],confs[i]["syn_out"])) for i in range(len(confs))]	# Inversion Profiles 1
	else:	
		fits = [p.read_profile(os.path.join(confs[i]["path"],confs[i]["inv_out"] + d.end_stokes)) for i in range(len(confs))]	# Inversion Profiles 1

	

	# Cut wave
	if confs[0]['mode'] == "1C" or confs[0]["mode"] == "2C":
		[obss[i].cut_to_wave(confs[i]["range_wave"]) for i in range(len(confs))]
	elif confs[0]['mode'] == "MC":
		num = fits[0].indx[0]
		if "-num" in sys.argv:
			num = int(sys.argv[sys.argv.index("-num")+1])
		ind1 = 0

		for n in range(len(confs)):
			for i in range(len(confs[n]["atoms"])):
				if str(num) in confs[n]["atoms"][i]:
					ind1 = i

					# define map as it is used later for saving
					confs[n]["map"] = [0,0,0,0]
					obss[n].cut_to_wave(np.array([confs[n]["range_wave"][ind1]]))
					fits[n].cut_to_wave(np.array([confs[n]["range_wave"][ind1]]))
	# Write in absolute wavelengths
	elif confs[0]["mode"] == "SY":
		for i in range(len(confs)):
			fits[i].transform_wave_sir_to_abs(os.path.join(confs[i]["path"],confs[i]["lines"]))
			fits[i].data_cut_wave = True
			confs[i]["map"] = [0,0,0,0]

	# Observation from synthesis
	if confs[0]["mode"] != "SY":
		lls, Is, Qs, Us, Vs = [obss[i].wave for i in range(len(confs))], [obss[i].stki[x[i],y[i]] for i in range(len(confs))],[obss[i].stkq[x[i],y[i]] for i in range(len(confs))],[obss[i].stku[x[i],y[i]] for i in range(len(confs))],[obss[i].stkv[x[i],y[i]]  for i in range(len(confs))]
	
	# Load fit data
	lls_fit, I_fits, Q_fits, U_fits, V_fits = [fits[i].wave for i in range(len(confs))], [fits[i].stki[x[i],y[i]] for i in range(len(confs))],[fits[i].stkq[x[i],y[i]] for i in range(len(confs))],[fits[i].stku[x[i],y[i]] for i in range(len(confs))],[fits[i].stkv[x[i],y[i]] for i in range(len(confs))]

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
	if confs[0]['mode'] == "MC":
		ll0 = sir.determine_line_core(os.path.join(confs[0]["path"],confs[0]["line"]),num)

	label_x = 0

	temp = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if temp != '0':
		if confs[0]["mode"] != "SY":
			lls = [lls[i]-float(temp) for i in range(len(lls))]
		else:
			lls = [lls_fit[i]-float(temp) for i in range(len(lls_fit))]
		lls_fit = [lls_fit[i]-float(temp) for i in range(len(lls_fit))]
		label_x = temp
	else:
		# Keep the same wavelengths as in the Grid file
		label_x = str(ll0)
		if confs[0]["mode"] != "SY":
			lls = [lls[i]-ll0 for i in range(len(lls))]
		else:
			lls = [lls_fit[i]-ll0 for i in range(len(lls_fit))]
		lls_fit = [lls_fit[i]-ll0 for i in range(len(lls_fit))]


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
	if confs[0]["mode"] == "MC":
		llabel = "Syn. Profile"
	else:
		llabel = "Obs. Profile"
	if confs[0]["mode"] != "SY":
		ax1.plot(lls[0], Is[0], "x", label=llabel)
		ax2.plot(lls[0], Qs[0], "x", label=llabel)
		ax3.plot(lls[0], Us[0], "x", label=llabel)
		ax4.plot(lls[0], Vs[0], "x", label=llabel)

	if confs[0]["mode"] != "SY":
		for i in range(len(confs)):
			ax1.plot(lls[i], I_fits[i], "-", label = "Best Fit (" + labels[i] + ")")
			ax2.plot(lls[i], Q_fits[i], "-", label = "Best Fit (" + labels[i] + ")")
			ax3.plot(lls[i], U_fits[i], "-", label = "Best Fit (" + labels[i] + ")")
			ax4.plot(lls[i], V_fits[i], "-", label = "Best Fit (" + labels[i] + ")")
	else:
		for i in range(len(confs)):
			ax1.plot(lls_fit[i], I_fits[i], "-", label = labels[i])
			ax2.plot(lls_fit[i], Q_fits[i], "-", label = labels[i])
			ax3.plot(lls_fit[i], U_fits[i], "-", label = labels[i])
			ax4.plot(lls_fit[i], V_fits[i], "-", label = labels[i])

	# Set xlimits
	ax1.set_xlim(lls[0][0], lls[0][-1])
	ax2.set_xlim(lls[0][0], lls[0][-1])
	ax3.set_xlim(lls[0][0], lls[0][-1])
	ax4.set_xlim(lls[0][0], lls[0][-1])

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
		ax1.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots
		plt.tight_layout(pad=2,h_pad=0.0)
	else:
		ax2.legend(bbox_to_anchor=(1.01,0.95))
		# set the spacing between subplots	
		plt.tight_layout(pad=2.5)

	plt.savefig(savepath + "inversion_multiple_stokes_x" + str(x) + "_y" + str(y) + add)

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

	for j in range(len(inputs)):
		if inputs[j] in sys.argv:
			fig, ax1 = plt.subplots(figsize=(8,7))

			# Add label
			ax1.scatter([], [], color="w", alpha=0, label=add_label)

			# '0C5DA5', 'FF2C00', '00B945', 'FF9500', '845B97', '474747', 'C20078','054907'
			# Plot
			# Synthesised model = Real model
			if confs[0]["mode"] == "MC" or confs[0]['mode'] == "2C":
				if confs[0]['mode'] == "2C":
					llabel = "Best Fit M. ("
					for i in range(len(confs)):
						ax1.plot(syns[i].tau, syns[i].get_attribute(inputs[m][1:])[x[i],y[i]], label=f"{llabel}{labels[i]})")
				else:
					llabel = "Syn Model"
					ax1.plot(syns[0].tau, syns[0].get_attribute(inputs[m][1:])[x[i],y[i]], label=f"{llabel}")
					ax1.plot(syns[0].tau, syns[0].get_attribute(inputs[m][1:])[x[i],y[i]], label=f"{llabel}")
	
			if confs[0]['mode'] == "2C":
					llabel = "Best Fit M. 2 ("
			else:
					llabel = "Best Fit ("
			for i in range(len(confs)):
				ax1.plot(fits[i].tau, fits[i].get_attribute(inputs[m][1:])[x[i],y[i]], label=f"{llabel}{labels[i]})")

			# Set xlimits
			ax1.set_xlim(fits[0].tau[0], fits[0].tau[-1])

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{c}$")
			ax1.set_ylabel(labels[index[j]])

			if index[j] == 2 or index[j] == 9:
				ax1.semilogy()

			# Legend
			ax1.legend()
		
			ax1.set_title(titles[i])
			# set the spacing between subplots
			plt.tight_layout(pad=2)
			plt.savefig(savepath + "inversion_multiple_x" + str(x) + "_y"  + str(y) + str(inputs[j][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = syns[0].tau[-1]
	if confs[0]["mode"] == "MC" or confs[0]['mode'] == "2C":
		[syns[i].set_limit(lim_max) for i in range(len(confs))]
		[syns[i].set_limit(lim_max) for i in range(len(confs))]
	[phys[i].set_limit(lim_max) for i in range(len(confs))]

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0))
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

	if confs[0]["mode"] == "MC" or confs[0]['mode'] == "2C":
		if confs[0]["mode"] == "2C":
			llabel = "Best Fit M. 1"
			for i in range(len(confs)):
				ax1.plot(syns[i].tau, syns[i].T[x[i],y[i]], label=f"{llabel} ({labels[i]})")
				ax2.plot(syns[i].tau, syns[i].B[x[i],y[i]], label=f"{llabel} ({labels[i]})")
				ax3.plot(syns[i].tau, syns[i].vlos[x[i],y[i]], label=f"{llabel} ({labels[i]})")
				ax4.plot(syns[i].tau, syns[i].gamma[x[i],y[i]], label=f"{llabel} ({labels[i]})")
		else:
			llabel = "Syn. Model"
			ax1.plot(syns[0].tau, syns[0].T[x[i],y[i]], label=f"{llabel}")
			ax2.plot(syns[0].tau, syns[0].B[x[i],y[i]], label=f"{llabel}")
			ax3.plot(syns[0].tau, syns[0].vlos[x[i],y[i]], label=f"{llabel}")
			ax4.plot(syns[0].tau, syns[0].gamma[x[i],y[i]], label=f"{llabel}")

	if confs[0]['mode'] == "2C":
		llabel = "Best Fit M. 2"
	else:
		llabel = "Best Fit"
	for i in range(len(confs)):
		ax1.plot(phys[i].tau, phys[i].T[x[i],y[i]], label=f"{llabel} ({labels[i]})")	
		ax2.plot(phys[i].tau, phys[i].B[x[i],y[i]], label=f"{llabel} ({labels[i]})")
		ax3.plot(phys[i].tau, phys[i].vlos[x[i],y[i]], label=f"{llabel} ({labels[i]})")
		ax4.plot(phys[i].tau, phys[i].gamma[x[i],y[i]], label=f"{llabel} ({labels[i]})")

	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(phys[0].tau[0], lim_max)
	ax2.set_xlim(phys[0].tau[0], lim_max)
	ax3.set_xlim(phys[0].tau[0], lim_max)
	ax4.set_xlim(phys[0].tau[0], lim_max)

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

	  
	plt.savefig(savepath + "inversion_multiple_result_x" + str(x) + "_y"  + str(y) + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	# Config files
	if "," in sys.argv[1]:
		confs = sys.argv[1].replace(", ", ",").split(',')
		confs = [sir.read_config(i, check=False) for i in confs]
		temp = sys.argv[1].split(",")
		for i in range(len(temp)):
			# Correct for ./ paths
			if confs[i]['path'][:2] == "./":
				confs[i]['path'] = temp[i][:temp[i].rfind('/')+1] + confs[i]['path'][2:]
	else:
		confs  = [sir.read_config(sys.argv[1])]
	
	if "," in sys.argv[2]:
		x = sys.argv[2].replace(", ", ",").split(',')
		x = [int(i) for i in x]
	else:
		x  = [int(sys.argv[2])]

	if "," in sys.argv[3]:
		y = sys.argv[3].replace(", ", ",").split(',')
		y = [int(i) for i in y]
	else:
		y  = [int(sys.argv[3])]

	#Labels
	labels = []
	if ("-labels" in sys.argv):
		labels = sys.argv[sys.argv.index("-labels")+1].split(",")
	else:
		for i in range(len(confs)):
			labels.append(input(f"Label {i+1}: "))
	
	inversion_multiple(confs, x, y, labels)




