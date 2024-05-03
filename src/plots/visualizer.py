#!/usr/bin/env python3
"""
Visualizes the spectra and model interactively at different pixels.
"""
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(sys.path[0] + "/..")
from matplotlib.backend_bases import MouseButton
import sir
import obs
from os.path import exists
import definitions as d
import signal
import profile_stk as p
import model as m

def signal_handling(signum,frame):
	"""
	Closes the program without printing the traceback as this is the way to finsish the program.
	"""
	print("\rClosing visualizer ...")
	if 'fig' in globals():
		plt.close(fig)
	if 'fig1' in globals():
		plt.close(fig1)
	if 'fig2' in globals():
		plt.close(fig2)
	sys.exit()


def help():
	print("visualizer - Visualizes the spectra and model.")
	print("Usage: python visualizer [OPTION]")
	print()
	sir.option("[1. Pos.]","Config")
	sir.option("[2. Pos.]","Wavelength at which stuff is plotted in A or log tau value if a phys. model is plotted as the map")
	print()
	sir.option("-data","Rel. path to the spectral veil corrected data if standard labelling is not used, optional.")	
	sir.option("-stokes","Rel. path to the Stokes result if standard labelling is not used, optional.")
	sir.option("-models","Rel. path to the Models of the inversion if standard labelling is not used.")
	sir.option("-chi","Rel. path to the chi2 file of the inversion if standard labelling is not used.")
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")
	sir.option("-label","Add label text (optional)")
	sir.option("-title","Title in Stokes plot")
	sir.option("-T","Plot temperature in K")
	sir.option("-Pe","Plot electron pressure in dyn/cm^2")
	sir.option("-vmicro","Plot microturbulence in cm/s")
	sir.option("-B","Plot magentic field strength in Gauss")
	sir.option("-vlos","Plot line of sight velocity in km/s")
	sir.option("-inc","Plot inclination in deg")
	sir.option("-azi","Plot azimuth in deg")
	sir.option("-z","Plot real height in km")
	sir.option("-Pg","Plot gas pressure in dyn/cm^2")
	sir.option("-rho","Plot density")
	sir.option("-vertical","Plot spectra vertically")
	sir.option("-chi2","Plot chi2")
	sir.option("-wave_V","Plot Stokes V in another wavelength position in A")

	sys.exit()

def check_range(wave_inv, wave):
	"""
	Check if the given wavelength is in the range

	Parameter
	---------
	wave_inv : numpy 1d array
		Array containing the wavelength from the inversion binary file
	wave : float
		Wavelength in A to be checked

	Return
	------
	Wavelength in the range
	Wavelength index in the range

	"""
	# Check for the wavelength if it is in the inversion range:
	if wave < np.min(wave_inv):
		print("Note that the given wavelength is not in the inversion range. Take closest one in range...")
		wave = np.min(wave_inv)
		return wave, 0
	elif wave > np.max(wave_inv):
		print("Note that the given wavelength is not in the inversion range. Take closest one in range...")
		wave = np.max(wave_inv)
		return wave, len(wave_inv)-1

	return wave, np.argmin(abs(wave_inv-wave))


def determine_resolution(string, pos, primary=False):
	"""
	Cuts a string from xrandr and extracts the width and height of the monitor

	Parameter
	---------
	string : str
		String where the information is extracted
	pos : int
		Position in the string where the information is extracted
	primary : bool, optional
		Determines wheter the chosen position is the primary monitor or not. Default: False
	"""
	if primary:
		width, height = string[pos][1:].split("+")[0].split()[1].split("x")
		added_width = int(string[pos][1:].split("+")[1])
		height = int(height)
	else:
		width, height = string[pos].split("+")[0].split('x')
		added_width = int(string[pos][1:].split("+")[1])
		height = int(height)
	return [added_width,int(width)], height


def get_display_size_2monitors():
	"""
	Gets the size of the display which is needed for positioning the figure for the second monitor

	Parameter
	---------
	None
	
	Return
	------
	width : int
		Width of the screen in pixel
	height : int
		Height of the screen in pixel
	"""
	import subprocess

	cmd = ['xrandr']
	cmd2 = ['grep', ' connected']
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	p2 = subprocess.Popen(cmd2, stdin=p.stdout, stdout=subprocess.PIPE)
	p.stdout.close()

	resolution_string, junk = p2.communicate()
	resolution = str(resolution_string).split("connected")
	
	width = height = width2 = height2 = -1
	# 3 Monitors
	if len(resolution) == 4:
		if "primary" in resolution[2]:
			width, height = determine_resolution(resolution, 2, True)
			width2, height2 = determine_resolution(resolution, 3, False)
		elif "primary" in resolution[3]:
			width, height = determine_resolution(resolution, 3, True)
			width2, height2 = determine_resolution(resolution, 2, False)
		else:
			width, height = determine_resolution(resolution, 1, True)
			width2, height2 = determine_resolution(resolution, 2, False)
	# 2 Monitors
	elif len(resolution) == 3:
		if "primary" in resolution[1]:
			width, height = determine_resolution(resolution, 1, True)
			width2, height2 = determine_resolution(resolution, 2, False)
		elif "primary" in resolution[2]:
			width, height = determine_resolution(resolution, 2, True)
			width2, height2 = determine_resolution(resolution, 1, False)
	# Do sth to at least assign existing values
	else:
		if "primary" not in resolution[1]:
			width, height = determine_resolution(resolution, 1, False)
		if "primary" not in resolution[2]:
			width2, height2 = determine_resolution(resolution, 2, False)
	
	return width, height, width2, height2


def get_display_size():
	"""
	Gets the size of the display which is needed for positioning the figure

	Parameter
	---------
	None
	
	Return
	------
	width : int
		Width of the screen in pixel
	height : int
		Height of the screen in pixel
	"""
	try:
		import subprocess
		cmd = ['xrandr']
		cmd2 = ['grep', '*']
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
		p2 = subprocess.Popen(cmd2, stdin=p.stdout, stdout=subprocess.PIPE)
		p.stdout.close()

		resolution_string, junk = p2.communicate()
		resolution = str(resolution_string).split()
		
		width, height = resolution[1].split('x')
		width = [int(width), int(width)]
		# Consider two or three monitors
		width2 = height2 = -1
		if len(resolution) > 5:
			width, height, width2, height2 = get_display_size_2monitors()
	except:
		print("The function 'get_display_size_2monitors' did not work properly.")
		print("Consider creating an issue on the gitlab page with your output of xrandr.")


	return width, int(height), width2, int(height2)

def move_figure(position, width = 10, height = 10, Monitor = "1"):
	'''
	Move and resize a window to a set of standard positions on the screen.
	Possible positions are:
	top, bottom, left, right, top-left, top-right, bottom-left, bottom-right
	To be used: define px (abs. width) and py (abs. height) as a global variable with the result from get_display_size()

	Source: https://stackoverflow.com/questions/7449585/how-do-you-set-the-absolute-position-of-figure-windows-with-matplotlib

	Parameter
	---------
	position : string
		Set the figure at this position.
	width : int, optional
		Width of the figure in pixel. Default: 10 (use it if one does not use right1 and right2)
	height : int, optional
		Height of the figure in pixel. Default: 10 (use it if one does not use right1 and right2)
	Monitor : string, optional
		Plot the figure on the second monitor
		
	Return
	------
	None	

	'''
	mgr = plt.get_current_fig_manager()
	
	aw = int(width)
	ah = int(height)
	
	
	# px etc. are defined with px, py, px2, py2 = get_display_size()
	if Monitor == "1":
		if position == "top":
			# x-top-left-corner, y-top-left-corner, x-width, y-width (in pixels)
			mgr.window.setGeometry(int(px[0]/2-aw/2),0, aw,ah)
		elif position == "bottom":
			mgr.window.setGeometry(int(px[0]/2-aw/2),int(py), aw,ah)
		elif position == "left":
			mgr.window.setGeometry(0, 0, int(px[1]/2),py)
		elif position == "right":
			mgr.window.setGeometry(int(px[0]), int(py/2+ah/2), aw,ah)
		elif position == "right1":
			mgr.window.setGeometry(int(px[0]), 0, int(px[1]/2),int(py/2.2))
		elif position == "right2":
			mgr.window.setGeometry(int(px[0]), int(py/2), int(px[1]/2),int(py/2.2))
		elif position == "top-left":
			mgr.window.setGeometry(0, 0, aw,ah)
		elif position == "top-right":
			mgr.window.setGeometry(px[0], 0, aw,ah)
		elif position == "bottom-left":
			mgr.window.setGeometry(0, py, aw,ah)
		elif position == "bottom-right":
			mgr.window.setGeometry(int(px[0]), int(py), aw,ah)

	elif Monitor == "2":
		if position == "top":
			# x-top-left-corner, y-top-left-corner, x-width, y-width (in pixels)
			mgr.window.setGeometry(int(px2[0]/2-aw/2),0, aw,ah)
		elif position == "bottom":
			mgr.window.setGeometry(int(px2[0]/2-aw/2),int(py2), aw,ah)
		elif position == "left":
			mgr.window.setGeometry(0, 0, int(px[1]/2),py2)
		elif position == "right":
			mgr.window.setGeometry(int(px2[0] - (aw/2+0.5)), int(py2/2+ah/2), aw,ah)
		elif position == "right1":
			mgr.window.setGeometry(int(px2[0]), 0, int(px2[1]),int(py2/2.2))
		elif position == "right2":
			mgr.window.setGeometry(int(px2[0]), int(py2/2), int(px2[1]),int(py2/2.2))
		elif position == "top-left":
			mgr.window.setGeometry(0, 0, aw,ah)
		elif position == "top-right":
			mgr.window.setGeometry(px2[0], 0, aw,ah)
		elif position == "bottom-left":
			mgr.window.setGeometry(0, py2, aw,ah)
		elif position == "bottom-right":
			mgr.window.setGeometry(int(px2[0]), int(py2), aw,ah)

def on_click_1C(event, obs, fit, Models, Map):
	"""
	Opens a figure with the I, Q, U and V in the chosen pixel in the picture.
	The chose pixel is determined by clicking on a figure

	Parameter
	---------
	event : unknown
 		Event containing the data of the chosen pixel
	obs : profile_stk
 		Observations as a profile_stk class
	fit : profile_stk
 		Inversion results as a profile_stk class
	Models : Model
 		Class with the model of the inversions
	Map : array
		Map of the image used in extent

	"""
	plt.ion() # No Error from QApplication
	if event.button is MouseButton.LEFT:
		xlims = False
		global fig1
		global fig2
		global ax1
		global ax2
		global ax3
		global ax4
		global ax11
		global ax21
		global ax31
		global ax41

		if (event.xdata == None):
			return
		elif (event.ydata == None):
			return
		
		# Read and print the position on the Map	
		x = int(event.xdata + 0.5) - Map[0]
		y = int(event.ydata + 0.5) - Map[2]
		print('[STATUS] Plot data at (%i,%i)' % (x + Map[0],y + Map[2]))


		# Set xlimits as before and also the size
		if 'fig1' in globals():
			xlims = True
			xlim1 = ax1.get_xlim()
			xlim2 = ax2.get_xlim()
			xlim3 = ax3.get_xlim()
			xlim4 = ax4.get_xlim()
			ax1.cla()
			ax2.cla()
			ax3.cla()
			ax4.cla()

		if 'fig2' in globals():
			ax11.cla()
			ax21.cla()
			ax31.cla()
			ax41.cla()

		# Set figure size as before
		if not xlims:
			fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
			plt.get_current_fig_manager().set_window_title('Stokes_' + str(rand))
			# In this position, the width and height are ignored which were set on the windows!
			size1 = fig1.get_size_inches()
			if px2 == -1:
				move_figure("right1")
			else:
				move_figure("right1", Monitor="2")

		ll1 = obs.wave
		ll2 = fit.wave

		ax1.plot(ll1, obs.stki[x,y,:], '-', label='Observation')
		ax1.plot(ll2, fit.stki[x,y,:], linestyle='dotted', label='Fit')

		ax2.plot(ll1, obs.stkq[x,y,:], '-')
		ax2.plot(ll2, fit.stkq[x,y,:], linestyle='dotted')

		ax3.plot(ll1, obs.stku[x,y,:], '-')
		ax3.plot(ll2, fit.stku[x,y,:], linestyle='dotted')

		ax4.plot(ll1, obs.stkv[x,y,:], '-')
		ax4.plot(ll2, fit.stkv[x,y,:], linestyle='dotted')

		if xlims:
			ax1.set_xlim(xlim1)
			ax2.set_xlim(xlim2)
			ax3.set_xlim(xlim3)
			ax4.set_xlim(xlim4)
		else:
			ax1.set_xlim(np.min(ll2), np.max(ll2))
			ax2.set_xlim(np.min(ll2), np.max(ll2))
			ax3.set_xlim(np.min(ll2), np.max(ll2))
			ax4.set_xlim(np.min(ll2), np.max(ll2))
		
		ax1.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax2.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax3.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax1.set_ylabel(r'$I/I_c$')
		ax2.set_ylabel(r'$Q/I_c$')
		ax3.set_ylabel(r'$U/I_c$')
		ax4.set_ylabel(r'$V/I_c$')

		fig1.suptitle(f"Pixel ({x},{y})")
		ax1.legend(fontsize=12)
		fig1.tight_layout()
		
		if not xlims:
			fig2, ((ax11,ax21),(ax31,ax41)) = plt.subplots(2,2)
			plt.get_current_fig_manager().set_window_title('Model_' + str(rand))
			size2 = fig2.get_size_inches()
			if px2 == -1:
				move_figure("right2")
			else:
				move_figure("right2", Monitor="2")
			


		log_tau = Models.log_tau
		ind_max = np.argmin(np.abs(log_tau	+ 3))
		
		log_tau = log_tau[:ind_max+1]
		ax11.plot(log_tau,Models.T[x,y,:ind_max+1])
		ax21.plot(log_tau,Models.B[x,y,:ind_max+1])
		ax31.plot(log_tau,Models.vlos[x,y,:ind_max+1])
		ax41.plot(log_tau,Models.gamma[x,y,:ind_max+1])
		
		ax11.set_title("Temperature T")
		ax21.set_title("Magnetic Field B")
		ax31.set_title(r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$")
		ax41.set_title(r"Inclination $\gamma$")

		ax11.set_xlabel(r"$\log \tau$")
		ax21.set_xlabel(r"$\log \tau$")
		ax31.set_xlabel(r"$\log \tau$")
		ax41.set_xlabel(r"$\log \tau$")

		ax11.set_ylabel(r"T [K]")
		ax21.set_ylabel(r"B [G]")
		ax31.set_ylabel(r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$")
		ax41.set_ylabel(r"$\gamma$ [deg]")



		ax11.set_xlim(log_tau[0],log_tau[-1])
		ax21.set_xlim(log_tau[0],log_tau[-1])
		ax31.set_xlim(log_tau[0],log_tau[-1])
		ax41.set_xlim(log_tau[0],log_tau[-1])

		fig2.suptitle(f"Pixel ({x},{y})")

		fig2.tight_layout()
		
		plt.show()

def visualizer_1C(conf, wave):
	"""
	Plots the Stokes vector with obs. and fits and the model depending on which pixel is clicked on in the appearing figure.

	Parameter
	---------
	config : dict
		All the information from the config file
	wave : float
		Wavelength in A where the Stokes Map is plotted / Or tau when model is plotted in the map

	"""
	# Import library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/../mml.mplstyle'):
		plt.style.use(dirname + '/../mml.mplstyle')
	elif "mml" in plt.style.available:
		plt.style.use('mml')

	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	path = conf["path"]
	Map = conf['map']

	if "-data" not in sys.argv:
		stokes = p.read_profile(os.path.join(conf["path"], conf['cube_inv']))
		stokes.cut_to_map(conf["map"])
	else:
		filename = sys.argv[sys.argv.index("-data")+1]
		stokes = p.read_proile(filename)
		stokes.cut_to_map(conf["map"])

	if "-stokes" not in sys.argv:
		stokes_inv = p.read_profile(os.path.join(path,conf['inv_out'] + d.end_stokes))
	else:
		filename = sys.argv[sys.argv.index("-stokes")+1]
		stokes_inv = p.read_profile(filename)
	
	if "-models" not in sys.argv:
			models_inv = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models))
	else:
			filename = sys.argv[sys.argv.index("-models")+1]
			models_inv = m.read_model(filename)
	
	# Check whether a model or Ic should be plotted
	use_model = False
	if ("-T" in sys.argv or "-Pe" in sys.argv or "-vmicro" in sys.argv or
		"-B" in sys.argv or "-vlos" in sys.argv or "-inc"  in sys.argv or
		"-azi" in sys.argv or "-z" in sys.argv or "-Pg" in sys.argv or
		"-rho" in sys.argv or "-chi2" in sys.argv):
		use_model = True
	
	# wave is used as tau position
	if use_model:
		tau = wave
		wave = 0 # Print the first value in the wavelength range
	waves_inv = stokes_inv.wave
	

	print("Opening visualizer ...")
	global px, py, px2, py2
	px, py, px2, py2 = get_display_size()
	
	global fig

	# Random number between 0 and 99 so that plotted figures can be linked to each other if many are plotted
	global rand
	rand = int(np.random.uniform(0,100))

	if use_model:
		if "-chi" not in sys.argv:
			chi2_inv = np.load(os.path.join(path,conf['chi2']))
		else:
			filename = sys.argv[sys.argv.index("-chi")+1]
			chi2_inv = np.load(filename)

		ind = np.argmin(abs(models_inv.log_tau - tau))
		

		###################################################
		#			Plot physical parameters			#
		###################################################

		# Define labels and get from arguments which parameter should be plot
		inputs = ["___","-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-chi2", "-Bz"]
		labels = ["", "T [K]", r"$\log P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$\log P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$",r"$\chi^2$", r"$B \cdot \cos \gamma$ [G]"]
		titles   = ["",r"Temperature T", r"Electron Pressure $\log P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength B",
					r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
				r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $\log P_g$", r"Density $\rho$", r"$\chi^2$", r"Line-of-Sight Magnetic Field $B \cdot \cos \gamma$"]
		cmap = [None,None,None,None,None,'seismic','jet','hsv',None,None,None,'gist_gray', None]
		limits = [[None,None],[np.min(models_inv[:,:,1,ind]),np.max(models_inv[:,:,1,ind])],[None,None],[None,None],[None,None],[-5,5],[0,180],[0,180],[None, None],[None, None],[None, None],[None,None],[-2000,2000]]
		i = 0
		imgscale = 7
		frac = (Map[3]+1-Map[2])/(Map[1]+1-Map[0])
		for i in range(len(inputs)):
			if inputs[i] in sys.argv:
				# Plot
				fig, ax = plt.subplots(figsize=[imgscale * 1.4, imgscale * frac])
				plt.get_current_fig_manager().set_window_title('Image_' + str(rand))
				size = fig.get_size_inches()
				move_figure("left", int(size[0]*fig.dpi),int(size[1]*fig.dpi))
				ax.set_title(titles[i] + r" @ $\log \tau = $" + str(tau))
				if inputs[i] == '-chi2':
					ax.set_title(titles[i])
					im = ax.imshow(chi2_inv.transpose(), cmap=cmap[i], origin = d.origin, vmin = limits[i][0], vmax = limits[i][1], extent=Map)	
				elif inputs[i] == '-Bz':
					ax.set_title(titles[i])
					im = ax.imshow((models_inv.B*np.cos(models_inv.gamma/180*np.pi)).transpose(), cmap=cmap[i], origin = d.origin, vmin = limits[i][0], vmax = limits[i][1], extent=Map)					
				else:
					im = ax.imshow(models_inv.get_attribute(inputs[i][2:])[:,:,ind].transpose(), cmap=cmap[i], origin = d.origin, vmin = limits[i][0], vmax = limits[i][1], extent=Map)
				
				# Set labels
				ax.set_xlabel(r"x [Pixels]")
				ax.set_ylabel(r"y [Pixels]")

				############
				# Colorbar #
				cbar = fig.colorbar(im, ax=ax, fraction=0.046 * frac, pad=0.04)
				cbar.set_label(label = labels[i], loc = 'center')
				############
				break
		plt.connect('button_press_event', lambda event: on_click_1C(event, stokes, stokes_inv, models_inv, Map))

		plt.show()

	else:
		# Check whether the wavelength is in range
		wave, wave_ind = check_range(stokes_inv.wave, wave)
		fig, ax = plt.subplots()
		plt.get_current_fig_manager().set_window_title('Image_' + str(rand))
		size = fig.get_size_inches()
		move_figure("left", int(size[0]*fig.dpi),int(size[1]*fig.dpi))

		im = ax.imshow(stokes_inv.stki[:,:,wave_ind].transpose(), #cmap = 'gist_gray',
					origin=d.origin, extent=Map)
		plt.connect('button_press_event', lambda event: on_click_1C(event, stokes, stokes_inv, models_inv, Map))
		plt.show()


def on_click_2C(event, obs, fit, Models1, Models2, Map):
	"""
	Opens a figure with the I, Q, U and V in the chosen pixel in the picture.
	The chose pixel is determined by clicking on a figure

	Parameter
	---------
	event : unknown
 		Event containing the data of the chosen pixel
	obs : profile_stk
 		Observations as a profile_stk class
	fit : profile_stk
 		Inversion results as a profile_stk class
	Models1 : Model
 		Class with the 1st models of the inversions
	Models2 : Model
 		Class with the 2nd models of the inversions
	Map : array
		Map of the image used in extent

	"""
	plt.ion() # No Error from QApplication
	if event.button is MouseButton.LEFT:
		xlims = False
		global fig1
		global fig2
		global ax1
		global ax2
		global ax3
		global ax4
		global ax11
		global ax21
		global ax31
		global ax41
		
		if event.xdata == None:
			return
		
		# Read and print the position on the Map	
		x = int(event.xdata + 0.5) - Map[0]
		y = int(event.ydata + 0.5) - Map[2]
		print('[STATUS] Plot data at (%i,%i)' % (x + Map[0],y + Map[2]))


		# Set xlimits as before and also the size
		if 'fig1' in globals():
			xlims = True
			xlim1 = ax1.get_xlim()
			xlim2 = ax2.get_xlim()
			xlim3 = ax3.get_xlim()
			xlim4 = ax4.get_xlim()
			ax1.cla()
			ax2.cla()
			ax3.cla()
			ax4.cla()

		if 'fig2' in globals():
			ax11.cla()
			ax21.cla()
			ax31.cla()
			ax41.cla()

		# Set figure size as before
		if not xlims:
			fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
			plt.get_current_fig_manager().set_window_title('Stokes_' + str(rand))
			if px1 == -1:
				move_figure("right1")
			else:
				move_figure("right1", Monitor="2")

		ll1 = obs.wave
		ll2 = fit.wave
		ax1.plot(ll1, obs.stki[x,y,:], '-', label="Observation")
		ax1.plot(ll2, fit.stki[x,y,:], linestyle='dotted', label="Fit")

		ax2.plot(ll1, obs.stkq[x,y,:], '-')
		ax2.plot(ll2, fit.stkq[x,y,:], linestyle='dotted')

		ax3.plot(ll1, obs.stku[x,y,:], '-')
		ax3.plot(ll2, fit.stku[x,y,:], linestyle='dotted')

		ax4.plot(ll1, obs.stkv[x,y,:], '-')
		ax4.plot(ll2, fit.stkv[x,y,:], linestyle='dotted')

		if xlims:
			ax1.set_xlim(xlim1)
			ax2.set_xlim(xlim2)
			ax3.set_xlim(xlim3)
			ax4.set_xlim(xlim4)
		else:
			ax1.set_xlim(np.min(ll2), np.max(ll2))
			ax2.set_xlim(np.min(ll2), np.max(ll2))
			ax3.set_xlim(np.min(ll2), np.max(ll2))
			ax4.set_xlim(np.min(ll2), np.max(ll2))
		
		ax1.legend(fontsize=12)

		ax1.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax2.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax3.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax4.set_xlabel(r'$\Delta \lambda$ [\AA]', loc='center')
		ax1.set_ylabel(r'$I/I_c$')
		ax2.set_ylabel(r'$Q/I_c$')
		ax3.set_ylabel(r'$U/I_c$')
		ax4.set_ylabel(r'$V/I_c$')
		fig1.suptitle(f"Pixel ({x},{y})")
		plt.tight_layout()

		if not xlims:		
			fig2, ((ax11,ax21),(ax31,ax41)) = plt.subplots(2,2)
			plt.get_current_fig_manager().set_window_title('Model_' + str(rand))
			if px2 == -1:
				move_figure("right2")
			else:
				move_figure("right2", Monitor="2")

		log_tau = Models1.log_tau
		ind_max = np.argmin(np.abs(log_tau	+ 3))
		
		log_tau = log_tau[:ind_max+1]
		ax11.plot(log_tau,Models1.T[x,y,:ind_max+1], label = "Model 1")
		ax21.plot(log_tau,Models1.B[x,y,:ind_max+1])
		ax31.plot(log_tau,Models1.vlos[x,y,:ind_max+1])
		ax41.plot(log_tau,Models1.gamma[x,y,:ind_max+1])
		
		ax11.plot(log_tau,Models2.T[x,y,:ind_max+1], label = "Model 2")
		ax21.plot(log_tau,Models2.B[x,y,:ind_max+1])
		ax31.plot(log_tau,Models2.vlos[x,y,:ind_max+1])
		ax41.plot(log_tau,Models2.gamma[x,y,:ind_max+1])

		ax11.legend()

		ax11.set_title("Temperature T")
		ax21.set_title("Magnetic Field B")
		ax31.set_title(r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$")
		ax41.set_title(r"Inclination $\gamma$")

		ax11.set_xlabel(r"$\log \tau$")
		ax21.set_xlabel(r"$\log \tau$")
		ax31.set_xlabel(r"$\log \tau$")
		ax41.set_xlabel(r"$\log \tau$")

		ax11.set_ylabel(r"T [K]")
		ax21.set_ylabel(r"B [G]")
		ax31.set_ylabel(r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$")
		ax41.set_ylabel(r"$\gamma$ [deg]")



		ax11.set_xlim(log_tau[0],log_tau[-1])
		ax21.set_xlim(log_tau[0],log_tau[-1])
		ax31.set_xlim(log_tau[0],log_tau[-1])
		ax41.set_xlim(log_tau[0],log_tau[-1])
		fig2.suptitle(f"Pixel ({x},{y})")

		plt.tight_layout()
		
		plt.show()

def visualizer_2C(conf, wave):
	"""
	Plots the Stokes vector with obs. and fits and the model depending on which pixel is clicked on in the appearing figure.

	Parameter
	---------
	conf : dict
		All the information from the config file
	wave : float
		Wavelength in A where the Stokes Map is plotted / Or tau when model is plotted in the map

	"""
	# Import library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/../mml.mplstyle'):
		plt.style.use(dirname + '/../mml.mplstyle')
	elif "mml" in plt.style.available:
		plt.style.use('mml')

	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	path = conf["path"]
	Map = conf['map']

	# Check whether a model or Ic should be plotted
	use_model = False
	if ("-T" in sys.argv or "-Pe" in sys.argv or "-vmicro" in sys.argv or
		"-B" in sys.argv or "-vlos" in sys.argv or "-inc"  in sys.argv or
		"-azi" in sys.argv or "-z" in sys.argv or "-Pg" in sys.argv or
		"-rho" in sys.argv or "-chi2" in sys.argv):
		use_model = True
	
	# wave is used as tau position
	if use_model:
		tau = wave
		wave = 0 # Print the first value in the wavelength range

	if "-data" not in sys.argv:
		stokes = p.read_profile(os.path.join(conf["path"], conf['cube_inv']))
		stokes.cut_to_map(conf["map"])
	else:
		filename = sys.argv[sys.argv.index("-data")+1]
		stokes = p.read_proile(filename)
		stokes.cut_to_map(conf["map"])

	if "-stokes" not in sys.argv:
		stokes_inv = p.read_profile(os.path.join(path,conf['inv_out'] + d.end_stokes))
	else:
		filename = sys.argv[sys.argv.index("-stokes")+1]
		stokes_inv = p.read_profile(filename)
	
	if "-models1" not in sys.argv:
			models1_inv = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models1))
	else:
			filename = sys.argv[sys.argv.index("-models")+1]
			models1_inv = m.read_model(filename)
	if "-models2" not in sys.argv:
			models2_inv = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models2))
	else:
			filename = sys.argv[sys.argv.index("-models")+1]
			models2_inv = m.read_model(filename)

	waves = stokes_inv.wave
	# Check whether the wavelength is in range
	wave, wave_ind = check_range(waves, wave)

	print("Opening visualizer ...")
	global px, py, px2, py2
	px, py, px2, py2 = get_display_size()

	global fig

	# Random number between 0 and 99 so that plotted figures can be linked to each other if many are plotted
	global rand
	rand = int(np.random.uniform(0,100))
	

	if use_model:
		if "-chi" not in sys.argv:
			chi2_inv = np.load(os.path.join(path,conf['chi2']))
		else:
			filename = sys.argv[sys.argv.index("-chi")+1]
			chi2_inv = np.load(filename)

		ind = np.argmin(abs(models1_inv[0,0,0,:] - tau))
			
		###################################################
		#			Plot physical parameters			#
		###################################################

		# Define labels and get from arguments which parameter should be plot
		inputs = ["___","-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-chi2", "-Bz"]
		labels = ["", "T [K]", r"$\log P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", "z [km]", r"$\log P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$",r"$\chi^2$", r"$B \cdot \cos \gamma$ [G]"]
		titles   = ["",r"Temperature T", r"Electron Pressure $\log P_e$",r"Microturbulence Velocity $\mathrm{v}_{\mathrm{micro}}$", r"Magnetic Field Strength B",
					r"Line-of-Sight Velocity $\mathrm{v}_{\mathrm{los}}$",
				r"Inclination $\gamma$", r"Azimuth $\phi$", r"Height $z$", r"Gas Pressure $\log P_g$", r"Density $\rho$", r"$\chi^2$", r"Line-of-Sight Magnetic Field $B \cdot \cos \gamma$"]
		cmap = [None,None,None,None,None,'seismic','jet','hsv',None,None,None,'gist_gray', None]
		limits = [[None,None],[np.min(models1_inv[:,:,1,ind]),np.max(models1_inv[:,:,1,ind])],[None,None],[None,None],[0,4500],[-5,5],[0,180],[0,180],[None, None],[None, None],[None, None],[None,None],[-2000,2000]]
		i = 0
		imgscale = 7
		frac = (Map[3]+1-Map[2])/(Map[1]+1-Map[0])
		for i in range(len(inputs)):
			if inputs[i] in sys.argv:
				# Plot
				fig, ax = plt.subplots(figsize=[imgscale * 1.4, imgscale * frac])
				plt.get_current_fig_manager().set_window_title('Image_' + str(rand))
				size = fig.get_size_inches()
				move_figure("left", int(size[0]*fig.dpi),int(size[1]*fig.dpi))
				ax.set_title(titles[i] + r" @ $\log \tau = $" + str(tau))
				if inputs[i] == '-chi2':
					ax.set_title(titles[i])
					im = ax.imshow(chi2_inv.transpose(), cmap=cmap[i], origin = d.origin, vmin = limits[i][0], vmax = limits[i][1], extent=Map)
				if inputs[i] == '-Bz':
					ax.set_title(titles[i])
					im = ax.imshow((models1_inv.B*np.cos(models1_inv.gamma/180*np.pi)).transpose(), cmap=cmap[i], origin = d.origin, vmin = limits[i][0], vmax = limits[i][1], extent=Map)					
				else:
					im = ax.imshow(models1_inv.get_attribute(inputs[i][2:])[:,:,ind].transpose(), cmap=cmap[i], origin = d.origin, vmin = limits[i][0], vmax = limits[i][1], extent=Map)
					
				# Set labels
				ax.set_xlabel(r"x [Pixels]")
				ax.set_ylabel(r"y [Pixels]")

				############
				# Colorbar #
				cbar = fig.colorbar(im, ax=ax, fraction=0.046 * frac, pad=0.04)
				cbar.set_label(label = labels[i], loc = 'center')
				############
				break
		plt.connect('button_press_event', lambda event: on_click_2C(event, stokes, stokes_inv, models1_inv, models2_inv, Map))

		plt.show()

	else:
		fig, ax = plt.subplots()
		plt.get_current_fig_manager().set_window_title('Image_' + str(rand))
		size = fig.get_size_inches()
		move_figure("left", int(size[0]*fig.dpi),int(size[1]*fig.dpi))

		im = ax.imshow(stokes_inv.stki[:,:,wave_ind].transpose(), #cmap = 'gist_gray',
					origin=d.origin, extent=Map)
		plt.connect('button_press_event', lambda event: on_click_2C(event, stokes, stokes_inv, models1_inv, models2_inv, Map))

		plt.show()

# Used if executed directly
if __name__ == "__main__":
	signal.signal(signal.SIGINT,signal_handling) 
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])
	if conf["mode"] == "1C":
		visualizer_1C(conf, float(sys.argv[2]))
	elif conf["mode"] == "2C":
		visualizer_2C(conf, float(sys.argv[2]))
	elif conf["mode"] == "MC":
		print("[visualizer] Not defined for mode 'MC'")
	else:
		print("[visualizer] Unknown mode")