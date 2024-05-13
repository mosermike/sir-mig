"""
Plots the chi2 values of multiple chi2 values
"""

import numpy as np 
import sys, os
import matplotlib.pyplot as plt
from os.path import exists
sys.path.append(os.path.join(os.path.dirname(__file__), "../..")) 
import sir
import definitions as d


def help():
	"""
	Help Page
	"""
	print("plot_chi2 - Plots the chi2 values in a histogram")
	print("Usage: python plot_chi2 [OPTION]")
	print()
	sir.option("[1. Pos]","Configs, given as an python array")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
	sir.option("-title","Title in Histogram")
	sir.option("-max","Biggest bin number")
	
	sys.exit()


def chi2(confs, labels):
	"""
	Plot the chi2 histogram for different simulations
	
	Parameter
	---------
	confs : list
		Dictionaries of the config file
	labels : list
		String list with the labels corresponding to the list above
	"""
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

	# Read Chi2 files
	chis = []
	for i in range(len(confs)):
		chis.append(np.load(os.path.join(confs[i]['path'], confs[i]['chi2'])))

	# Define stuff for the plots from input
	if "-max" in sys.argv:
		bins = np.linspace(0, float(sys.argv[sys.argv.index("-max")+1]), 100)
	else:
		bins = np.linspace(0, d.chi2_lim, 100)

	title = ''
	if "-title" in sys.argv:
		title = sys.argv[sys.argv.index("-title")+1]
	label = '_'
	if "-label" in sys.argv:
		label = sys.argv[sys.argv.index("-label")+1]

	add = ''
	if "-add" in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	save = ''
	if "-save" in sys.argv:
		save = sys.argv[sys.argv.index("-save")+1]

	# Plot the histogram
	fig, ax = plt.subplots()

	ax.scatter([], [], color="w", alpha=0, label=label)
	for i in range(len(chis)):
		ax.hist(chis[i], histtype='step', label=labels[i], bins=bins)
	ax.set_xlabel(r"$\chi^2$")
	ax.set_ylabel(r"Entries")

	ax.set_xlim(0, bins[-1])

	if title != '':
		ax.set_title(title)

	ax.legend()

	# set the spacing between subplots
	plt.tight_layout()

	plt.savefig(save + "chi2" + add)


# Used if executed directly
if __name__ == "__main__":
	############################
	# READ INPUT AND LOAD DATA #
	############################
	if "-h" in sys.argv:
		help()

	# Config files
	if "," in sys.argv[1]:
		confs = sys.argv[1].replace(", ", ",").strip('][').split(',')
		confs = [sir.read_config(i, check=False) for i in confs]
	else:
		confs = [sir.read_config(sys.argv[1])]
	# Labels
	labels = []
	for i in range(len(confs)):
		ll = input(f"Label {i+1}: ")
		labels.append(ll)

	chi2(confs, labels)
