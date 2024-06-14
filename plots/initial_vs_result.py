import sys
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(sys[0] +  "/../src") 
import sir
import definitions as d
from os.path import exists


def help():
	"""
	Help Page
	"""
	print("plot_initial_vs_result - Plots the initial value of the guess model and the result at a to-be-specified log tau value")
	print("Usage: python plot_initial_vs_result.py [OPTION]")
	print()
	sir.option("[1. Pos]","List of Guess files")
	sir.option("[2. Pos]","List of result models files")
	sir.option("[3. Pos]","List of chi2 files")
	sir.option("[4. Pos]","At what log tau value it is evaluated")

	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
	sir.option("-title","Title in Stokes plot")
	
	sys.exit()

def plot_initial_vs_result():
	"""
	Plots the initial value vs the resulting value for different physical parameters and the chi2 value in a 2x2 plot

	Parameters
	----------
	None

	Returns
	-------
	None

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


	##########################
	#	Define variables	#
	##########################
	Guess_file  = sys.argv[1].split(",") # List of the Guess Models
	Result_file = sys.argv[2].split(",") # List of the Results Models
	inv_file = sys.argv[3].split(",") # List of the chi2 files

	arg = float(sys.argv[4])	# at what log_tau it is evaluated

	# Additional labels
	save = ''
	if "-save" in sys.argv:
		save = sys.argv[sys.argv.index("-save")+1]
		# Create savepath
		if not exists(save[0:save.rfind("/")]):
			os.mkdir(save[0:save.rfind('/')])

	# Additional text in saved figures
	add = ''
	if "-add" in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Title of the last plot
	title = ''
	if "-title" in sys.argv:
		title = sys.argv[sys.argv.index("-title")+1]

	# Title x position
	xtitle = 0.5
	if '-xtitle' in sys.argv:
		xtitle = float(sys.argv[sys.argv.index("-xtitle")+1])

	# Log tau must be everywhere the same for simplicity
	log_tau, T, _,_,B,vlos,inc,_,_,_,_ = sir.read_model(Guess_file[0]) # Read first file for shapes

	guess  = np.zeros(shape=(len(Guess_file),5,len(log_tau))) # Parameters from the guesses
	result = np.copy(guess) # Result Models
	error  = np.copy(guess) # Errors from the best fit
	chi2   = np.copy(guess) # Chi2 from the best fit
	Num = len(Guess_file)
	for i in range(Num):
		# Guess
		log_tau,T,_,_,B,vlos,inc,_,_,_,_ = sir.read_model(f"{Guess_file[i]}")
		guess[i][0] = log_tau
		guess[i][1] = T
		guess[i][2] = B
		guess[i][3] = vlos/1e5
		guess[i][4] = inc
		# Result
		log_tau,T,_,_,B,vlos,inc,_,_,_,_ = sir.read_model(f"{Result_file[i]}")
		result[i][0] = log_tau
		result[i][1] = T
		result[i][2] = B
		result[i][3] = vlos/1e5
		result[i][4] = inc
		# Error
		log_tau,T,_,_,B,vlos,inc,_,_,_,_ = sir.read_model(f"{Result_file[i]}")
		error[i][0] = log_tau
		error[i][1] = T
		error[i][2] = B
		error[i][3] = vlos/1e5
		error[i][4] = inc
		chi2[i] = sir.read_chi2(inv_file[i])
		
	# Determine index, where log_tau = arg
	ind = np.where(log_tau == arg)[0][0]

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16*0.8,12*0.8), sharex=True)

	put_labels = True
	n = 1


	# Hollow circle for the first list
	ax1.plot(np.arange(Num)+1, guess[:,1, ind], 'o', markerfacecolor='none', markeredgecolor='#FF2C00', markersize=8, label="Initial Guess", color='#FF2C00')
	# Filled circle for the second list
	ax1.plot(np.arange(Num)+1, result[:,1, ind], 'o', markerfacecolor='#0C5DA5', markersize=8, label="Result", color = '#0C5DA5')

	# Hollow circle for the first list
	ax2.plot(np.arange(Num)+1, guess[:,2, ind], 'o', markerfacecolor='none', markeredgecolor='#FF2C00', markersize=8, label="Initial Guess", color='#FF2C00')
	# Filled circle for the second list
	ax2.plot(np.arange(Num)+1, result[:,2, ind], 'o', markerfacecolor='#0C5DA5', markersize=8, label="Result", color = '#0C5DA5')

	# Hollow circle for the first list
	ax3.plot(np.arange(Num)+1, guess[:,3, ind], 'o', markerfacecolor='none', markeredgecolor='#FF2C00', markersize=8, label="Initial Guess", color='#FF2C00')
	# Filled circle for the second list
	ax3.plot(np.arange(Num)+1, result[:,3, ind], 'o', markerfacecolor='#0C5DA5', markersize=8, label="Result", color = '#0C5DA5')

	ax4.plot(np.arange(Num)+1, chi2, "o")



	ax4.set_xlabel(r"Repetition", loc = "center")
	ax3.set_xlabel(r"Repetition", loc = "center")

	titles   = [r"Temperature T", r"Magnetic Field Strength B",
		       r"Line-of-Sight Velocity $v_{\mathrm{los}}$", r"$\chi^2$"]

	labels = [r"T ($\log \tau_5 = " + str(arg) + r"$)[K]",r"B ($\log \tau_5 = " + str(arg) + r"$)[G]", r"$v_{\mathrm{los}}$ ($\log \tau_5 = " + str(arg) + r"$) $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", r"$\chi^2$"]

	ax1.set_ylabel(labels[0])
	ax2.set_ylabel(labels[1])
	ax3.set_ylabel(labels[2])
	ax4.set_ylabel(labels[3])

	ax1.legend(frameon=True, loc="best")

	#ax1.set_title(titles[0])
	#ax2.set_title(titles[3])
	#ax3.set_title(titles[4])
	#ax4.set_title(titles[5])

	if title != '':
		fig.suptitle(title, y=0.98, x=xtitle)
	plt.tight_layout(pad=1.5)
	plt.savefig(save + 'init_vs_res' + add)

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	plot_initial_vs_result()

