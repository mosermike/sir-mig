"""
Plots the chi2 values of multiple chi2 values
"""

import numpy as np 
import sys, os
import matplotlib.pyplot as plt
from os.path import exists
sys.path.append(os.path.join(os.path.dirname(__file__), "../../src")) 
import sir
import definitions as d
import chi2_stk as c

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


def chi2(confs : list, labels : list):
	"""
	Plot the chi2 histogram for different simulations
	
	Parameter
	---------
	confs : list
		Dictionaries of the config file
	labels : list
		String list with the labels corresponding to the list above

	Return
	------
	None

	Other Parameters
	----------------

	-max [float]
		Maximum chi2 bin
	-title [str]
		Title of the plot
	-save [str]
		Path where the figure is saved
	-add [str]
		Additional string at the end of the figure file name
	"""
	# Import library
	sir.mpl_library()

	# Read Chi2 files
	chis = []
	for i in range(len(confs)):
		chis.append(c.read_chi2(os.path.join(confs[i]['path'], confs[i]['chi2'])))

	# Define stuff for the plots from input
	if "-max" in sys.argv:
		bins = np.linspace(0, float(sys.argv[sys.argv.index("-max")+1]), 100)
	else:
		bins = np.linspace(0, d.chi2_lim, 100)

	title = ''
	if "-title" in sys.argv:
		title = sys.argv[sys.argv.index("-title")+1]

	add = ''
	if "-add" in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	save = ''
	if "-save" in sys.argv:
		save = sys.argv[sys.argv.index("-save")+1]

	# Plot the histogram
	fig, ax = plt.subplots()

	#ax.scatter([], [], color="w", alpha=0, label=label)
	for i in range(len(chis)):
		ax.hist(chis.tot[i], histtype='step', label=labels[i], bins=bins)
	ax.set_xlabel(r"$\chi_{\text{tot}}^2$")
	ax.set_ylabel(r"Entries")

	ax.set_xlim(0, bins[-1])

	if title != '':
		ax.set_title(title)

	ax.legend()

	# set the spacing between subplots
	plt.tight_layout()

	plt.savefig(save + "chi2_tot" + add)

	# Plot the histogram
	fig, ax = plt.subplots()

	#ax.scatter([], [], color="w", alpha=0, label=label)
	for i in range(len(chis)):
		ax.hist(chis.stki[i], histtype='step', label=labels[i], bins=bins)
	ax.set_xlabel(r"$\chi_{I}^2$")
	ax.set_ylabel(r"Entries")

	ax.set_xlim(0, bins[-1])

	if title != '':
		ax.set_title(title)

	ax.legend()

	# set the spacing between subplots
	plt.tight_layout()

	plt.savefig(save + "chi2_stki" + add)

	# Plot the histogram
	fig, ax = plt.subplots()

	#ax.scatter([], [], color="w", alpha=0, label=label)
	for i in range(len(chis)):
		ax.hist(chis.stkq[i], histtype='step', label=labels[i], bins=bins)
	ax.set_xlabel(r"$\chi_Q^2$")
	ax.set_ylabel(r"Entries")

	ax.set_xlim(0, bins[-1])

	if title != '':
		ax.set_title(title)

	ax.legend()

	# set the spacing between subplots
	plt.tight_layout()

	plt.savefig(save + "chi2_stkq" + add)

	# Plot the histogram
	fig, ax = plt.subplots()

	#ax.scatter([], [], color="w", alpha=0, label=label)
	for i in range(len(chis)):
		ax.hist(chis.stku[i], histtype='step', label=labels[i], bins=bins)
	ax.set_xlabel(r"$\chi_U^2$")
	ax.set_ylabel(r"Entries")

	ax.set_xlim(0, bins[-1])

	if title != '':
		ax.set_title(title)

	ax.legend()

	# set the spacing between subplots
	plt.tight_layout()

	plt.savefig(save + "chi2_stku" + add)

	# Plot the histogram
	fig, ax = plt.subplots()

	#ax.scatter([], [], color="w", alpha=0, label=label)
	for i in range(len(chis)):
		ax.hist(chis.stkv[i], histtype='step', label=labels[i], bins=bins)
	ax.set_xlabel(r"$\chi^2$")
	ax.set_ylabel(r"Entries")

	ax.set_xlim(0, bins[-1])

	if title != '':
		ax.set_title(title)

	ax.legend()

	# set the spacing between subplots
	plt.tight_layout()

	plt.savefig(save + "chi2_stkv" + add)


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
		confs = [sir.read_config(i) for i in confs]
	else:
		confs = [sir.read_config(sys.argv[1])]
		
	# Labels
	labels = []
	for i in range(len(confs)):
		ll = input(f"Label {i+1}: ")
		labels.append(ll)

	chi2(confs, labels)
