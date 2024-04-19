import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import sir
import definitions as d

def help():
	"""
	Print the help page about how to use this script
	"""

	print("response - Plots a response function for maximal two lines in the Grid file.")
	print("Usage: python response.py [OPTION]")
	print()
	sir.option("[1. Pos]","Filename of the response function text file from SIR")
	sir.option("[2. Pos]","Grid file")
	sir.option("[3. Pos]","Line file")
	sir.option("[4. Pos]","The Actual Model")
	print()
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
	
	sys.exit()

#############################################################################

def load_data(filename):
	"""
	Loads the data from the response function file

	Parameter
	---------
	filename : str
		String containing the filename of the response function

	Return
	------
	ntau : numbers of tau
	nlam : numbers of lambda
	I : Stokes I
	Q : Stokes Q
	U : Stokes U
	V : Stokes V

	"""
	# Open the file in read mode
	with open(filename, 'r') as file:
		lines = file.readlines()

	# Extracting variables from the file
	temp = list(filter(None, lines[0].replace('\n','').split(' '))) # Remove newlines and spaces
	ntau, nlam4 = map(int, temp)

	# Calculate n
	n = ntau * nlam4

	# Initialize data array
	data = np.zeros(n)

	# Read data from the file
	data_line = lines[1].replace("  "," ").replace('\n','').split(' ')
	data_line = list(filter(None, data_line)) # Remove empty strings
	data = np.array(list(map(float, data_line)))

	# Reshape the data array
	nlam = int(nlam4 // 4)
	d2 = data.reshape((ntau, 4, nlam))

	# Extracting arrays
	I = d2[:, 0, :]
	Q = d2[:, 1, :]
	U = d2[:, 2, :]
	V = d2[:, 3, :]
	return ntau, nlam, I, Q, U, V

#############################################################################

def break_two_axes(ax1,ax2, frac):
	"""
	Breaks to axis and sets the limits (tested for relation 1:2, 2:1 and 1:1

	Parameter
	---------
	ax1 : matplotlib axis
		First axis
	ax2 : matplotlib axis
		Second axis
	frac : float
		Fraction defining the relation between axis 1 and axis 2

	Return
	------
	None

	"""
	# Hide the spines between ax1 and ax2
	ax1.spines['right'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	
	# No ticks between ax1 and ax2
	ax1.yaxis.tick_left()
	ax2.set_yticks([])
    
	d = .02  # how big to make the diagonal lines in axes coordinates

	# Put the diagonal lines according to the fraction (may be adjusted if it is wrong)
	if frac < 1:
		kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
		ax2.plot((-4.5*d, 0.5*d), (1-d, 1+d), **kwargs)
		ax2.plot((-2*d, +3*d), (1-d, 1+d), **kwargs)
		ax2.plot((-4.5*d, 0.5*d), (-d, +d), **kwargs)
		ax2.plot((-2*d, +3*d), (-d, +d), **kwargs)
	elif frac > 1:
		kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
		ax2.plot((-2*d, 0), (1-d, 1+d), **kwargs)
		ax2.plot((-d, +d), (1-d, 1+d), **kwargs)

		#kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
		ax2.plot((-2*d, 0), (-d, d), **kwargs)
		ax2.plot((-d, +d), (-d, d), **kwargs)
	else:
		kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
		ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
		ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

		kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
		ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
		ax2.plot((-d, +d), (-d, +d), **kwargs)

#############################################################################

def plot_response(filename, grid, line, modelname):
	"""
	Plots the response function


	Parameter
	---------
	filename : string
		Filename of the response function file
	grid : string
		Filename of the grid file
	line : string
		Filename of the line file
	modelname : string
		Filename of the used model
	"""

	# Import matplotlib library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/mml.mplstyle'):
		plt.style.use(dirname + '/mml.mplstyle')
	elif "mml" in plt.style.available:
		plt.style.use('mml')

	grid = sir.read_grid(grid)
	line = sir.read_line(line)

	# Read the grid
	Line      = grid['Line']
	Line_min  = grid['min']
	Line_step = grid['step']
	Line_max  = grid['max']
	points = []
	range_wave = []
	for i in range(len(Line)):	
		points.append(round((Line_max[i]+Line_step[i] - Line_min[i])/Line_step[i]+0.5))
		ll0 = line['wavelength'][np.where(line['Line'] == int(Line[i][0]))]
		range_wave.append([ll0[0]+Line_min[i]/1000,ll0[0]+Line_max[i]/1000])



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


	log_tau, T, Pe, vmicro, B, vlos, inc, azi, z, pg, rho = sir.read_model(modelname)
	ntau, nlam, I, Q, U, V = load_data(filename)

	##############################################
	#  Plot I, Q, U and V  at wave for obs		#
	##############################################
	# Change to relative wavelength to make the axis better readable
	rel = input('Put a wavelength to print in rel. wavelengths (0 = in absolute): ')

	range_wave = np.array(range_wave) - float(rel)

	# Defining minima and maxima for the colorbar and normalise
	I = I/np.max(I)
	I_min = np.min(I)
	Q = Q/np.max(np.abs(Q))
	U = U/np.max(np.abs(U))
	V = V/np.max(np.abs(V))

	if len(points) == 2:
		Map1 = [range_wave[0,0],range_wave[0,1],log_tau[-1],log_tau[0]]
		Map2 = [range_wave[1,0],range_wave[1,1],nlam,log_tau[-1],log_tau[0]]

		frac = (range_wave[1,1]-range_wave[1,0])/(range_wave[0,1]-range_wave[0,0])

		imgscale = 9
		fig, axes = plt.subplots(2,5,figsize=(16,12),
		               gridspec_kw=dict(width_ratios=[4,4 * frac,2*(1+frac/5),4,4 * frac], wspace=0.025,hspace=0.5))


		# plot the same data on both axes
		im11 = axes[0,0].imshow(I[:,0:points[0]], vmin=I_min,vmax=1, cmap = 'Reds', extent=Map1, aspect='auto')
		im12 = axes[0,1].imshow(I[:,points[0]:] , vmin=I_min,vmax=1, cmap = 'Reds', extent=Map2, aspect='auto')
		im21 = axes[0,3].imshow(Q[:,0:points[0]], vmin=-1,vmax=1, cmap = 'RdBu', extent=Map1, aspect='auto')
		im22 = axes[0,4].imshow(Q[:,points[0]:] , vmin=-1,vmax=1, cmap = 'RdBu', extent=Map2, aspect='auto')
		im31 = axes[1,0].imshow(U[:,0:points[0]], vmin=-1,vmax=1, cmap = 'RdBu', extent=Map1, aspect='auto')
		im32 = axes[1,1].imshow(U[:,points[0]:] , vmin=-1,vmax=1, cmap = 'RdBu', extent=Map2, aspect='auto')
		im41 = axes[1,3].imshow(V[:,0:points[0]], vmin=-1,vmax=1, cmap = 'RdBu', extent=Map1, aspect='auto')
		im42 = axes[1,4].imshow(V[:,points[0]:] , vmin=-1,vmax=1, cmap = 'RdBu', extent=Map2, aspect='auto')


		axes[0,2].remove()
		axes[1,2].remove()

		break_two_axes(axes[0,0],axes[0,1], frac)
		break_two_axes(axes[0,3],axes[0,4], frac)
		break_two_axes(axes[1,0],axes[1,1], frac)
		break_two_axes(axes[1,3],axes[1,4], frac)

		#im1 = ax1.imshow(I/np.max(I), cmap = 'Reds', extent=Map, aspect='auto')
		#im2 = ax2.imshow(Q/np.max(np.abs(Q)), cmap = 'RdBu', extent=Map, aspect='auto')
		#im3 = ax3.imshow(U/np.max(np.abs(U)), cmap = 'RdBu', extent=Map, aspect='auto')
		#im4 = ax4.imshow(V/np.max(np.abs(V)), cmap = 'RdBu', extent=Map, aspect='auto')

		#####################
		#	Set labels	#
		#####################
		# Get to which parameter the RF are computed
		temp = ''
		if filename[-1] == 't':
			temp = r'$T$'
		elif filename[-1] == 'h':
			temp = r'$B$'
		elif filename[-3:] == 'inc':
			temp = r'$\gamma$'
		elif filename[-2:] == 'vz':
			temp = r'v$_{\mathrm{los}}$'
		if frac <= 1:
			axes[0,0].set_title(r'RF of Stokes $I$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
			axes[0,3].set_title(r'RF of Stokes $Q$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
			axes[1,0].set_title(r'RF of Stokes $U$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
			axes[1,3].set_title(r'RF of Stokes $V$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
		if frac > 1:
			axes[0,1].set_title(r'RF of Stokes $I$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
			axes[0,4].set_title(r'RF of Stokes $Q$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
			axes[1,1].set_title(r'RF of Stokes $U$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)
			axes[1,4].set_title(r'RF of Stokes $V$ to ' + f'{temp}', x = (4)/(4+4*frac), y = 1.01)

		#####################
		#	Set labels	#
		#####################
		rel_str = '' # String to be add to axis labels to symbolise the rel. wavelength
		if int(rel) != 0:
			rel_str = f"- {rel} "
		if "-vertical" not in sys.argv:
			axes[0,1].set_xlabel(r"$\lambda$ " + rel_str + "$[\AA]$")
			axes[0,4].set_xlabel(r"$\lambda$ " + rel_str + "$[\AA]$")
			axes[1,1].set_xlabel(r"$\lambda$ " + rel_str + "$[\AA]$")
		axes[1,4].set_xlabel(r"$\lambda$  " + rel_str + "$[\AA]$")

		axes[0,0].set_ylabel(r"$\log \tau_{500}$")
		axes[0,3].set_ylabel(r"$\log \tau_{500}$")
		axes[1,0].set_ylabel(r"$\log \tau_{500}$")
		axes[1,3].set_ylabel(r"$\log \tau_{500}$")

		frac = 2/frac
		############
		# Colorbar #
		cbar1 = fig.colorbar(im12, ax=axes[0,1], fraction=0.046 * frac, pad=0.04)
		cbar1.set_label(label = r'$I / I_c $', loc = 'center')
		cbar2 = fig.colorbar(im22, ax=axes[0,4], fraction=0.046 * frac, pad=0.04)
		cbar2.set_label(label = r'$Q / I_c $', loc = 'center')
		cbar3 = fig.colorbar(im32, ax=axes[1,1], fraction=0.046 * frac, pad=0.04)
		cbar3.set_label(label = r'$U / I_c $', loc = 'center')
		cbar4 = fig.colorbar(im42, ax=axes[1,4], fraction=0.046 * frac, pad=0.04)
		cbar4.set_label(label = r'$V / I_c $', loc = 'center')
		############

	# Only one line in Grid file
	elif len(points) == 1:
		Map = [range_wave[0,0],range_wave[0,1],log_tau[-1],log_tau[0]]

		frac = 1
		imgscale = 9
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=[imgscale * 1.4, imgscale * frac],
											#sharex = True, sharey = True
								)

		# plot the same data on both axes
		im1 = ax1.imshow(I, vmin=I_min,vmax=1, cmap = 'Reds', extent=Map, aspect='auto')
		im2 = ax2.imshow(Q, vmin=-1   ,vmax=1, cmap = 'RdBu', extent=Map, aspect='auto')
		im3 = ax3.imshow(U, vmin=-1   ,vmax=1, cmap = 'RdBu', extent=Map, aspect='auto')
		im4 = ax4.imshow(V, vmin=-1   ,vmax=1, cmap = 'RdBu', extent=Map, aspect='auto')

		#####################
		#	Set labels	#
		#####################
		temp = ''
		if filename[-1] == 't':
			temp = r'$T$'
		elif filename[-1] == 'h':
			temp = r'$B$'
		elif filename[-3:] == 'inc':
			temp = r'$\gamma$'
		elif filename[-2:] == 'vz':
			temp = r'v$_{\mathrm{los}}$'
		ax1.set_title(r'RF of Stokes $I$ to ' + f'{temp}')
		ax2.set_title(r'RF of Stokes $Q$ to ' + f'{temp}')
		ax3.set_title(r'RF of Stokes $U$ to ' + f'{temp}')
		ax4.set_title(r'RF of Stokes $V$ to ' + f'{temp}')

		#####################
		#	Set labels	#
		#####################
		rel_str = '' # String to be add to axis labels to symbolise the rel. wavelength
		if int(rel) != 0:
			rel_str = f"- {rel} "
		if "-vertical" not in sys.argv:
			ax1.set_xlabel(r"$\lambda$ " + rel_str + "[\AA]")
			ax2.set_xlabel(r"$\lambda$ " + rel_str + "[\AA]")
			ax3.set_xlabel(r"$\lambda$ " + rel_str + "[\AA]")
		ax4.set_xlabel(r"$\lambda$  " + rel_str + "[\AA]")

		ax1.set_ylabel(r"$\log \tau_{500}$")
		ax2.set_ylabel(r"$\log \tau_{500}$")
		ax3.set_ylabel(r"$\log \tau_{500}$")
		ax4.set_ylabel(r"$\log \tau_{500}$")

		############
		# Colorbar #
		cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.046 * frac, pad=0.04)
		cbar1.set_label(label = r'$I / I_c $', loc = 'center')
		cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046 * frac, pad=0.04)
		cbar2.set_label(label = r'$Q / I_c $', loc = 'center')
		cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.046 * frac, pad=0.04)
		cbar3.set_label(label = r'$U / I_c $', loc = 'center')
		cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.046 * frac, pad=0.04)
		cbar4.set_label(label = r'$V / I_c $', loc = 'center')
		############
	
	else:
		print('[NOTE] There are more than 2 lines in the Grid file. The x axis is not plotted with the real wavelength (not implemented)')

		Map = [0,len(I),log_tau[-1],log_tau[0]]

		frac = 1
		imgscale = 9
		fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=[imgscale * 1.4, imgscale * frac],
											#sharex = True, sharey = True
								)

		# plot the same data on both axes
		im1 = ax1.imshow(I, vmin=I_min,vmax=1, cmap = 'Reds', extent=Map, aspect='auto')
		im2 = ax2.imshow(Q, vmin=-1   ,vmax=1, cmap = 'RdBu', extent=Map, aspect='auto')
		im3 = ax3.imshow(U, vmin=-1   ,vmax=1, cmap = 'RdBu', extent=Map, aspect='auto')
		im4 = ax4.imshow(V, vmin=-1   ,vmax=1, cmap = 'RdBu', extent=Map, aspect='auto')

		#####################
		#	Set labels	#
		#####################
		temp = ''
		if filename[-1] == 't':
			temp = r'$T$'
		elif filename[-1] == 'h':
			temp = r'$B$'
		elif filename[-4:] == 'inc':
			temp = r'$\gamma$'
		elif filename[-3:] == 'vz':
			temp = r'v$_{\mathrm{los}}$'
		ax1.set_title(r'RF of Stokes $I$ to ' + f'{temp}')
		ax2.set_title(r'RF of Stokes $Q$ to ' + f'{temp}')
		ax3.set_title(r'RF of Stokes $U$ to ' + f'{temp}')
		ax4.set_title(r'RF of Stokes $V$ to ' + f'{temp}')

		#####################
		#	Set labels	#
		#####################
		if "-vertical" not in sys.argv:
			ax1.set_xlabel(r"Wavelength Point")
			ax2.set_xlabel(r"Wavelength Point")
			ax3.set_xlabel(r"Wavelength Point")
		ax4.set_xlabel(r"Wavelength Point")

		ax1.set_ylabel(r"$\log \tau_{500}$")
		ax2.set_ylabel(r"$\log \tau_{500}$")
		ax3.set_ylabel(r"$\log \tau_{500}$")
		ax4.set_ylabel(r"$\log \tau_{500}$")

		############
		# Colorbar #
		cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.046 * frac, pad=0.04)
		cbar1.set_label(label = r'$I / I_c $', loc = 'center')
		cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046 * frac, pad=0.04)
		cbar2.set_label(label = r'$Q / I_c $', loc = 'center')
		cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.046 * frac, pad=0.04)
		cbar3.set_label(label = r'$U / I_c $', loc = 'center')
		cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.046 * frac, pad=0.04)
		cbar4.set_label(label = r'$V / I_c $', loc = 'center')
		############
	

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title != "-1": # -1 => Print no title
		if "-vertical" in sys.argv:
			xtitle1 = 0.41
		else:
			xtitle1 = 0.5
		if title != '':
			fig.suptitle(title, y=0.98, x=xtitle1)

	#########################
	# Set Legend and Limits #
	#########################
	
	if len(points) == 1:
		plt.tight_layout(pad=2.0)

	temp = ''
	if filename[-1] == 't':
		temp = r'_T'
	elif filename[-1] == 'h':
		temp = r'_B'
	elif filename[-3:] == 'inc':
		temp = r'_inc'
	elif filename[-2:] == 'vz':
		temp = r'_vlos'
	plt.savefig(savepath + "response" + temp + add)
	
	return

#############################################################################

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	plot_response(sys.argv[1], sys.argv[2],sys.argv[3], sys.argv[4])



