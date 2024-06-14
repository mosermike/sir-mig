"""
Plots the result of the SIR inversion
"""

import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import sys
import os
sys.path.append(sys.path[0] + "/../src")
import os, sys, sir
import definitions as d
import profile_stk as p

def _help():
	print("Ic_test - Plots the continuum in the obs vs the fit to test the inversion.")
	print("Usage: python Ic_test.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Config")
	print()
	sir.option("-data","Rel. path to the spectral veil corrected data if standard labelling is not used, optional.")	
	sir.option("-stokes","Rel. path to the Stokes result if standard labelling is not used, optional.")
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-title","Title in plot")
	sir.option("-num","determines which range is used (1 = 1st line in grid, etc.) if not the first one is wanted.")

	sys.exit()

def Ic_test(conf, num):
	"""
	Plots Ic fit vs I obs in a plot. It uses the first element in the wavelength range as continuum

	Parameters
	----------
	config : dict
		Dict. with all the information from the config
	num : int
		determines which range is used (1 = 1st line in grid, etc.)
	
	Returns
	-------
	None

	Other Parameters
	----------------
	Additional parameters given as an argument when the script is executed.
	
	-data [str]
		Rel. path to the spectral veil corrected data if standard labelling is not used.
	-stokes [str]
		Rel. path to the Stokes result if standard labelling is not used.
	-save [str], optional
		Additional save path. Default './'.
	-add [str]
		Additional text in filenames.
	-title [str]
		Title in plot.
	-num [int]
		Determines which range is used (1 = 1st line in grid, etc.) if not the first one is wanted.
	"""
	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	path = conf["path"]
	stokes_inv = p.read_profile(os.path.join(path,conf['inv_out']))
	
	Map = conf['map']

	if "-data" not in sys.argv:
		stokes = p.read_profile(os.path.join(conf["path"],conf['cube']))
	else:
		filename = sys.argv[sys.argv.index("-data")+1]
		stokes = p.read_profile(filename)

	if "-data" not in sys.argv:
		stokes_inv = p.read_profile(os.path.join(path,conf['inv_out'] + d.end_stokes))
	else:
		filename = sys.argv[sys.argv.index("-stokes")+1]
		stokes_inv = p.read_profile(filename)
	
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


	# Determine the index to be plotted as continuum
	wave_ind1 = 0 # Index for inversion stokes
	wave_ind2 = np.argmin(np.abs(stokes.wave - stokes_inv.wave[0])) # Index for observations

	#############################################################
	#					PLOTTING STUFF					#
	#############################################################
	fig, ax = plt.subplots()

	############################
	# Plot the Stokes profiles #
	############################
	ax.plot(stokes.stki[Map[0]:Map[1]+1,Map[2]:Map[3]+1,wave_ind2].flatten(),stokes_inv.stki[:,:,wave_ind1].flatten(), '.', label="")
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	ax.plot([0,10],[0,10], '--', label="Expected Relation")
	ax.text(xlim[0]+0.05+0.01,xlim[0]+0.05, "Expected Relation",color='#FF2C00', rotation = 37.5)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)

	#####################
	#	Set labels	#
	#####################
	ax.set_xlabel(r"$\mathrm{I}_c^{\mathrm{obs}}$")
	ax.set_ylabel(r"$\mathrm{I}_c^{\mathrm{fit}}$")

	#####################
	#	Set title		#
	#####################

	if title != "-1":
		if title != '':
			ax.set_title(title)
		else:
			ax.set_title(r'Continuum Intensity $I_c$ Test')

	#########################
	# Set Legend and Limits #
	#########################
	plt.savefig(savepath + "inversion_test_ic" + add)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
	conf = sir.read_config(sys.argv[1])

	num = 0
	if '-num' in sys.argv:
		num = sys.argv[sys.argv.index("-num")+1]-1
	if conf['mode'] == "1C" or "2C":
		Ic_test(conf, int(num))
	else:
		print("[Ic_test] Mode unknown or not defined")




