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

def _help():
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

def synthesis_multiple(confs : list, x : list, labels : list):
	""""
	Plots the result of multiple inversions. The code assumes the same x and y position for all inversions. In addition, it also assumes that in mode "MC"
	the synthesis model is the same

	Parameters
	----------
	confs : list
		List with dictionaries of the config files
	x : list
		Model numbers (x positions)
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

	if(len(x) != len(confs)):
		raise RuntimeError("[inversion_multiple] Lengths of the arrays confs and x not the same!")
	
	if(np.all(np.array([i["mode"] for i in confs], dtype=str) != confs[0]["mode"])):
		raise ValueError(f"The mode must be the same for all configs. The first one is '{confs[0]['mode']}'.")
	
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################

	if confs[0]['mode'] == "MC":
		syns = [m.read_model(os.path.join(confs[i]["path"],confs[i]["syn_out"] + d.end_models)) for i in range(len(confs))]	# Synthesis Models 1
	elif confs[0]["mode"] == "SY":
		syns = [m.read_model(os.path.join(confs[i]["path"],confs[i]["syn_in"])) for i in range(len(confs))]
	
	if confs[0]['mode'] == "MC":
		obss = [p.read_profile(os.path.join(confs[i]["path"],confs[i]["syn_out"] + d.end_stokes)) for i in range(len(confs))]	# Synthesis Profiles 1
	else:
		obss = [p.read_profile(os.path.join(confs[i]["path"],confs[i]["syn_out"])) for i in range(len(confs))]

	

	# Cut wave
	if "-num" in sys.argv:
		num = int(sys.argv[sys.argv.index("-num")+1])
		ind1 = 0

		for n in range(len(confs)):
			confs[n]["map"] = [0,0,0,0]
			for i in range(len(confs[n]["atoms"])):
				if str(num) in confs[n]["atoms"][i]:
					ind1 = i
					obss[n].cut_to_wave(np.array([confs[n]["range_wave"][ind1]]))
		ll0 = sir.determine_line_core(os.path.join(confs[0]["path"],confs[0]["line"]),num)
	else:
		for i in range(len(confs)):
			obss[i].transform_wave_sir_to_abs(os.path.join(confs[i]["path"],confs[i]["line"]))
			obss[i].data_cut_wave = True
			confs[i]["map"] = [0,0,0,0]
	
	lls = [obss[i].wave for i in range(len(confs))]
	Is  = [obss[i].stki[x[i],0] for i in range(len(confs))]
	Qs  = [obss[i].stkq[x[i],0] for i in range(len(confs))]
	Us  = [obss[i].stku[x[i],0] for i in range(len(confs))]
	Vs  = [obss[i].stkv[x[i],0]  for i in range(len(confs))]
						
	
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
	

	label_x = 0

	temp = input("Put wavelength in A to which it should be relative (0 = change nothing): ")
	if temp != '0':
		lls = [lls[i]-float(temp) for i in range(len(lls))]
		label_x = temp
		
	linestyle_str = ['-', '--', ':', '-.', (0, (3, 5, 1, 5)), (0, (5, 10)),
              (0, (5, 10, 1, 10)), (0, (5, 2)), (0, (1, 1)), (0, (3, 5, 1, 5, 1, 5))]

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
	for i in range(len(confs)):
		ax1.plot(lls[i], Is[i], label=labels[i], linestyle=linestyle_str[i % len(linestyle_str)])
		ax2.plot(lls[i], Qs[i], label=labels[i], linestyle=linestyle_str[i % len(linestyle_str)])
		ax3.plot(lls[i], Us[i], label=labels[i], linestyle=linestyle_str[i % len(linestyle_str)])
		ax4.plot(lls[i], Vs[i], label=labels[i], linestyle=linestyle_str[i % len(linestyle_str)])

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
	
	if label_x != 0:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')#, loc='center')
		ax4.set_xlabel(r'$\Delta \lambda - $' + label_x + r' \AA')#, loc='center')
	else:
		if "-vertical" not in sys.argv:
			ax3.set_xlabel(r'$\lambda$ [\AA]')#, loc='center')
		ax4.set_xlabel(r'$\lambda$ [\AA]')#, loc='center')
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
		ax1.legend()
	else:
		ax1.legend()
	
	fig.savefig(savepath + "syn_multiple_stokes_x0_" + str(x[0]) + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["-T", '-Pe', '-vmicro', '-B', "-vlos", "-gamma", "-phi", "-z", "-Pg","-rho"]
	index  = [1,2,3,4,5,6,7,8,9,10]
	llabels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$"]
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
			for i in range(len(confs)):
				ax1.plot(syns[i].tau, syns[i].get_attribute(inputs[i][1:])[x[i],0], label=labels[i], linestyle = linestyle_str[i % len(linestyle_str)])
			
			# Set xlimits
			ax1.set_xlim(syns[0].tau[0], syns[0].tau[-1])

			# Set labels
			ax1.set_xlabel(r"$\log \tau_{c}$")
			ax1.set_ylabel(llabels[index[j]])

			if index[j] == 2 or index[j] == 9:
				ax1.semilogy()

			# Legend
			ax1.legend()
		
			ax1.set_title(titles[i])
			# set the spacing between subplots
			fig.savefig(savepath + "syn_multiple_x0_" + str(x[0]) + str(inputs[j][1:]) + add)
		
	# Plot T,B,vlos, inc in one figure
	lim_max = syns[0].tau[-1]
	[syns[i].set_limit(lim_max) for i in range(len(confs))]
	[syns[i].set_limit(lim_max) for i in range(len(confs))]
	

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(12,16), sharex=True,
			 gridspec_kw=dict(hspace=0), layout="compressed")
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True, layout="compressed")

	for i in range(len(confs)):
		ax1.plot(syns[i].tau, syns[i].T[x[i],0], label=labels[i], linestyle = linestyle_str[i % len(linestyle_str)])
		ax2.plot(syns[i].tau, syns[i].B[x[i],0], label=labels[i], linestyle = linestyle_str[i % len(linestyle_str)])
		ax3.plot(syns[i].tau, syns[i].vlos[x[i],0], label=labels[i], linestyle = linestyle_str[i % len(linestyle_str)])
		ax4.plot(syns[i].tau, syns[i].gamma[x[i],0], label=labels[i], linestyle = linestyle_str[i % len(linestyle_str)])
	
	#####################
	#	Set limits	#
	#####################	
	ax1.set_xlim(syns[0].tau[0], lim_max)
	ax2.set_xlim(syns[0].tau[0], lim_max)
	ax3.set_xlim(syns[0].tau[0], lim_max)
	ax4.set_xlim(syns[0].tau[0], lim_max)

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

	ax1.set_ylabel(llabels[1])
	ax2.set_ylabel(llabels[4])
	ax3.set_ylabel(llabels[5])
	ax4.set_ylabel(llabels[6])

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

	fig.savefig(savepath + "syn_multiple_result_x0_" + str(x[0]) + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
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
		nums = sys.argv[2].replace(", ", ",").split(',')
		nums = [int(i) for i in nums]
	else:
		nums  = [int(sys.argv[2])]
	#Labels
	labels = []
	if ("-labels" in sys.argv):
		labels = sys.argv[sys.argv.index("-labels")+1].split(",")
	else:
		for i in range(len(confs)):
			labels.append(input(f"Label {i+1}: "))
	
	synthesis_multiple(confs, nums, labels)




