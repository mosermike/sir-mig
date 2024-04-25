"""
Plots the result of the SIR inversion
"""

import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import os, sys
sys.path.append(sys.path[0] + "/../..")
sys.path.append(sys.path[0] + "/../../tools")
import sir, obs
import definitions as d
from change_config_path import change_config_path
import model as m
import profile_stk as p

# TODO do the figure size as in gris_firtez

def help():
	"""
	Help Page
	"""
	print("result - Plots the result of the inversion")
	print("Usage: python result.py [OPTION]")
	print()
	sir.option("1:","Config")
	sir.option("2:","wavelength at which stuff is plotted in A")
	sir.option("3:","tau value at which stuff is plotted")
	print()
	sir.option("-data:","Rel. path to the spectral veil corrected data if standard labelling is not used, optional.")	
	sir.option("-stokes:","Rel. path to the Stokes result if standard labelling is not used, optional.")
	sir.option("-models:","Rel. path to the Models of the inversion if standard labelling is not used.")
	sir.option("-errors:","Rel. path to the Errors of the inversion if standard labelling is not used.")
	sir.option("-chi:","Rel. path to the chi2 file of the inversion if standard labelling is not used.")
	sir.option("-save:","Additional save path (optional, default './')")
	sir.option("-add:","Additional text in filenames (optional)")
	sir.option("-label:","Add label text (optional)")
	sir.option("-title1:","Title in Result Stokes plot")
	sir.option("-title2:","Title in Obs. Stokes plot")
	sir.option("-title3:","Title in Model plot")
	sir.option("-xtitle:","Changing the x position of the title")
	sir.option("-T:","Plot temperature in K")
	sir.option("-Pe:","Plot electron pressure in dyn/cm^2")
	sir.option("-vmicro:","Plot microturbulence in cm/s")
	sir.option("-B:","Plot magentic field strength in Gauss")
	sir.option("-vlos:","Plot line of sight velocity in km/s")
	sir.option("-inc:","Plot inclination by subtracting in deg")
	sir.option("-azi:","Plot azimuth by adding in deg")
	sir.option("-z:","Plot real height in km")
	sir.option("-Pg:","Plot gas pressure in dyn/cm^2")
	sir.option("-rho:","Plot density")
	sir.option("-vertical:","Plot spectra vertically")
	sir.option("-chi2:","Plot chi2")
	sir.option("-plot_chi2:","Plot chi2 in the 4 subplots plot")
	sir.option("-waveV:","Plot Stokes V in another wavelength position in A")
	sir.option("-logT:","Plot Temperature at this log tau.")
	sir.option("-logB:","Plot Magnetic Field at this log tau.")
	sir.option("-logV:","Plot LoS Velocity at this log tau.")
	sir.option("-logI:","Plot Inclination at this log tau.")
	sir.option("-limitxy:","Limit in the x and y plot as a list xmin,xmax,ymin,xmax")
	sir.option("-limitT:","Set the limit for the colorbar in the temperature.")
	sir.option("-limitB:","Set the limit for the colorbar in the magnetic field.")
	sir.option("-limitchi2:","Set the limit for the colorbar in chi2.")
	sir.option("-limitI:","Set the limit for the colorbar in Stokes I.")
	sir.option("-arc:","Print x and y axis in arcseconds")
	sir.option("-flipx:","Mirror/Flip data as sometimes it is wrong in GRIS with the location on the sun")

	sys.exit()

def check_range(wave_inv, wave):
	"""
	Check if the given wavelength is in the range

	Parameter
	---------
	wave_inv : numpy array
		Array containing the wavelengths used in the inversion
	wave : float
		Wavelength in A to be checked

	Return
	------
	Wavelength in the range
	Index of the wavelength
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

	wave_ind = np.argmin(abs(wave_inv-wave))
	return wave, wave_ind


def plot(conf, wave, tau, waveV = -1):
	"""
	Plots the result of the inversion

	Parameter
	---------
	config : dict
		Dict. with all the information from the config
	wave : float
		Wavelength in A where the Stokes vector is plottet
	tau : float
		log tau value where the model is plotted
	waveV : float, optional
		Wavelength for Stokes V in A (-1 not used), Default: -1

	Return
	-------
	None

	"""

		# Import library
	dirname = os.path.split(os.path.abspath(__file__))[0]
	if exists(dirname + '/../../mml.mplstyle'):
		plt.style.use(dirname + '/../../mml.mplstyle')
	elif "mml" in plt.style.available:
		plt.style.use('mml')

	# Check if path exists
	if not exists(conf['path']):
		Inp = input("[NOTE] Path does not exist. You want to overwrite it with the actual path? [y/n] ")
		if Inp == "y":
			change_config_path(conf,os.path.abspath(os.getcwd()),)

	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	path = conf["path"]
	Map = conf['map']
		

	

	if "-data" not in sys.argv:
		stokes = p.read_profile(os.path.join(conf["path"],conf['cube_inv']))
	else:
		filename = sys.argv[sys.argv.index("-data")+1]
		stokes = p.read_profile(filename)

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
	if "-errors" not in sys.argv:
		errors_inv = m.read_model(os.path.join(path,conf['inv_out'] + d.end_errors))
	else:
		filename = sys.argv[sys.argv.index("-errors")+1]
		errors_inv = m.read_model(filename)
	
	if "-chi" not in sys.argv:
		chi2_inv = np.load(os.path.join(path,conf['chi2']))
	else:
		filename = sys.argv[sys.argv.index("-chi")+1]
		chi2_inv = np.load(filename)

	if "-logT" not in sys.argv:
		logT = tau
	else:
		logT = float(sys.argv[sys.argv.index("-logT")+1])

	if "-logB" not in sys.argv:
		logB = tau
	else:
		logB = float(sys.argv[sys.argv.index("-logB")+1])

	if "-logV" not in sys.argv:
		logV = tau
	else:
		logV = float(sys.argv[sys.argv.index("-logV")+1])

	if "-logI" not in sys.argv:
		logI = tau
	else:
		logI = float(sys.argv[sys.argv.index("-logI")+1])


	# Check whether the wavelength is in range
	wave, wave_ind = check_range(stokes_inv.wave, wave)
	
	if "-waveV" in sys.argv:
		waveV = float(sys.argv[sys.argv.index("-waveV")+1])
		waveV, waveV_ind = check_range(stokes_inv.wave, waveV)

	if int(waveV) != -1:
		waveV, waveV_ind = check_range(stokes_inv.wave, waveV)


	# Cut data in x and y position	
	if "-limitxy" in sys.argv:
		limit_xy = np.array([int(i) for i in sys.argv[sys.argv.index("-limitxy")+1].split(",")], dtype=int)

		# Cut data to the new range:
		#stokes = stokes[limit_xy[0]-Map[0]:limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2]:limit_xy[3]+1-(Map[3]+1)]
		stokes_inv = stokes_inv.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])
		models_inv = models_inv.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])
		errors_inv = errors_inv.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])

		# Save the new limits as the new Map
		Map = limit_xy

	waves = stokes.wave
	waves_inv = stokes_inv.wave

	if waveV == -1:
		waveV_ind1 = np.argmin(abs(waves-wave))
		waveV_ind2 = np.argmin(abs(waves_inv-wave))
	else:
		waveV_ind1 = np.argmin(abs(waves-waveV))
		waveV_ind2 = np.argmin(abs(waves_inv-waveV))
		waveV = waves[waveV_ind1]
	# Determine indexes for the wavelength
	wave_ind1 = np.argmin(abs(waves-wave))
	wave_ind2 = np.argmin(abs(waves_inv-wave))
	wave = waves[wave_ind1]

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = path + "/" + sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath[:savepath.rfind('/')]):
			os.mkdir(savepath[:savepath.rfind('/')])

	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Title
	title1 = ''
	if '-title1' in sys.argv:
		title1 = sys.argv[sys.argv.index("-title1")+1]
	title2 = ''
	if '-title2' in sys.argv:
		title2 = sys.argv[sys.argv.index("-title2")+1]
	title3 = ''
	if '-title3' in sys.argv:
		title3 = sys.argv[sys.argv.index("-title3")+1]

	# Additional label
	add_label = '_'
	if '-label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]

	# Restrict models to the given tau
	ind  = np.argmin(abs(models_inv.log_tau - tau))
	indT = np.argmin(abs(models_inv.log_tau - logT))
	indB = np.argmin(abs(models_inv.log_tau - logB))
	indV = np.argmin(abs(models_inv.log_tau - logV))
	indI = np.argmin(abs(models_inv.log_tau - logI))

	# Cut models to the specific log tau value
	models_inv.nval = ind
	models_inv.T = models_inv.T[:,:,:indT]
	models_inv.Pe = models_inv.Pe[:,:,:ind]
	models_inv.vmicro = models_inv.vmicro[:,:,:ind]
	models_inv.B = models_inv.B[:,:,:indB]
	models_inv.vlos = models_inv.vlos[:,:,:indV]
	models_inv.gamma = models_inv.gamma[:,:,:indI]
	models_inv.phi = models_inv.phi[:,:,:ind]
	models_inv.z = models_inv.z[:,:,:ind]
	models_inv.rho = models_inv.rho[:,:,:ind]
	models_inv.Pg = models_inv.Pg[:,:,:ind]


	taus = [tau for i in range(11)]
	taus[1] = logT
	taus[4] = logB
	taus[5] = logV
	taus[6] = logI
	
	# Pressure in log
	models_inv.Pe = np.log(models_inv.Pe)
	models_inv.Pg = np.log(models_inv.Pg)

	#############################################################
	#					PLOTTING STUFF					#
	#############################################################
	Markers = ["-", '--', 'dotted', 'dashdotdotted', 'densely dashed']

	linestyle_str = [
		'solid',	# Same as (0, ()) or '-'
		'dotted',    # Same as (0, (1, 1)) or ':'
		'dashed',    # Same as '--'
		'dashdot',
		'solid',	# Same as (0, ()) or '-'
		'dotted',    # Same as (0, (1, 1)) or ':'
		'dashed',    # Same as '--'
		'dashdot',
		(0, (3,10,1,10))
		] 


	stokes.cut_to_map(Map)
	I1 = stokes.stki
	Q1 = stokes.stkq
	U1 = stokes.stku
	V1 = stokes.stkv

	I2 = stokes_inv.stki
	Q2 = stokes_inv.stkq
	U2 = stokes_inv.stku
	V2 = stokes_inv.stkv

	limits_stokes1 = [  [np.min(I1[:,:,wave_ind1]), np.max(I1[:,:,wave_ind1])],
					[-np.max(np.abs(Q1[:,:,wave_ind1])), np.max(np.abs(Q1[:,:,wave_ind1]))],
					[-np.max(np.abs(U1[:,:,wave_ind1])), np.max(np.abs(U1[:,:,wave_ind1]))],
					[-np.max(np.abs(V1[:,:,wave_ind1])), np.max(np.abs(V1[:,:,waveV_ind1]))]
		]
	limits_stokes2 = [  [np.min(I2[:,:,wave_ind2]), np.max(I2[:,:,wave_ind2])],
					[-np.max(np.abs(Q2[:,:,wave_ind2])), np.max(np.abs(Q2[:,:,wave_ind2]))],
					[-np.max(np.abs(U2[:,:,wave_ind2])), np.max(np.abs(U2[:,:,wave_ind2]))],
					[-np.max(np.abs(V2[:,:,wave_ind2])), np.max(np.abs(V2[:,:,waveV_ind2]))]
		]


	if "-limitI" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitI")+1].split(',')
		limits_stokes1[0] = [int(i) for i in temp ]
		limits_stokes2[0] = [int(i) for i in temp ]

	##############################################
	#  Plot I, Q, U and V  at wave for obs		#
	##############################################
	if I2.shape[0] / I2.shape[1] >= 0.8:
		alpha = I2.shape[0] / 12
		figsize = [I2.shape[0] / alpha , I2.shape[1]/alpha-1]
		frac = figsize[1] / figsize[0]
	else:
		alpha = I2.shape[0] / 18
		figsize = [I2.shape[0] / alpha, I2.shape[1]/alpha-5]
		frac = figsize[1] / figsize[0] * 1.25
	
	####################################
	#	Plot settings for arcsecond	#
	####################################
	if "-arc" in sys.argv:
		infos = dict(np.loadtxt(os.path.join(path, d.header_infos), delimiter='=', dtype=str))
		if conf['instrument'] == 'GRIS':
			# y position starts at the end but I start with the pixels lower left
			y_pos = float(infos['CRVAL2']) + (stokes.shape[1]-1) * float(infos['CDELT2'])
			y_max = float(infos['CRVAL2']) # Used if flipx is used
			Map_plot = [float(infos['CRVAL1']) - float(infos['CDELT1']) * (Map[0]-1),
						float(infos['CRVAL1']) - float(infos['CDELT1']) * (Map[1]-1),
						y_pos - float(infos['CDELT2']) * (Map[2]-1),
						y_pos - float(infos['CDELT2']) * (Map[3]-1)
						]
		elif conf['instrument'] == 'Hinode':
			delta_y = float(infos['CDELT2'])  # Delta y of slit
			delta_x = float(infos['XSCALE'])
			x_pos = float(infos['XCEN'])
			y_pos = float(infos['YCEN'])
			y_max = y_pos + ((stokes.shape[1]-1) - infos['CRPIX2']) * delta_y # Used for flipx

			Map_plot = [
				(Map[0]-1) * delta_x + x_pos,
				(Map[1]-1) * delta_x + x_pos,
				(Map[2]-1 - infos['CRPIX2']) * delta_y + y_pos,
			    (Map[3]-1 - infos['CRPIX2']) * delta_y + y_pos
			]
	else:
		Map_plot = Map

	# Gris data is sometimes flipped along x! Perform real flip also in arcsec
	if "-flipx" in sys.argv:
		print("[NOTE]  Plots are flipped/mirrored along x")
		if d.origin == 'upper':
			origin = 'lower'
		else:
			origin = 'upper'
		Map_plot[2] = y_pos + (y_max - Map_plot[2])
		Map_plot[3] = y_pos + (y_max - Map_plot[3])
		Map_plot[2], Map_plot[3] = Map_plot[3], Map_plot[2]
	else:
		origin = d.origin

	if "-arc" in sys.argv:
		units = 'Arcsec'
	else:
		units = 'Pixels'


	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True,
			 gridspec_kw=dict(hspace=0), figsize=figsize)
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if (I2.shape[1] - I2.shape[0]) > 100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0],figsize[1]/4],
										layout="compressed",
									)
		else:
			fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=figsize,
										layout="compressed",
									)
		#fig.subplots_adjust(hspace=0.0)
	############################
	# Plot the Stokes profiles #
	############################
	im1 = ax1.imshow(I1[:,:,wave_ind1]  .transpose(), origin=origin, vmin = limits_stokes1[0][0], vmax = limits_stokes1[0][1], extent=Map_plot)
	im2 = ax2.imshow(Q1[:,:,wave_ind1]  .transpose(), origin=origin, vmin = limits_stokes1[1][0], vmax = limits_stokes1[1][1], cmap = 'PuOr', extent=Map_plot)
	im3 = ax3.imshow(U1[:,:,wave_ind1]  .transpose(), origin=origin, vmin = limits_stokes1[2][0], vmax = limits_stokes1[2][1], cmap = 'PuOr', extent=Map_plot)
	im4 = ax4.imshow(V1[:,:,waveV_ind1].transpose(), origin=origin, vmin = limits_stokes1[3][0], vmax = limits_stokes1[3][1], cmap = 'PuOr', extent=Map_plot)

	#####################
	#	Set labels	#
	#####################
	ax1.set_title(r'$\mathrm{I} / \mathrm{I}_c $ @' + "%.3f" % wave + r" \AA")
	ax2.set_title(r'$\mathrm{Q} / \mathrm{I}_c $')
	ax3.set_title(r'$\mathrm{U} / \mathrm{I}_c $')
	ax4.set_title(r'$\mathrm{V} / \mathrm{I}_c $')
	if waveV != -1:
		ax4.set_title(r'$\mathrm{V} / \mathrm{I}_c$ @' + "%.3f" % waveV + r" \AA")	


	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.057 * frac, pad=0.04)
	#cbar1.set_label(label = r'$I / I_c $', loc = 'center')
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.057 * frac, pad=0.04)
	#cbar2.set_label(label = r'$Q / I_c $', loc = 'center')
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.057 * frac, pad=0.04)
	#cbar3.set_label(label = r'$U / I_c $', loc = 'center')
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.057 * frac, pad=0.04)
	#cbar4.set_label(label = r'$V / I_c $', loc = 'center')
	############

	#####################
	#	Set labels	#
	#####################
	if "-vertical" not in sys.argv:
		ax3.set_xlabel(f"x [{units}]")
	ax4.set_xlabel(f"x [{units}]")
	ax1.set_ylabel(f"y [{units}]")
	ax3.set_ylabel(f"y [{units}]")

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title2 != "-1": # -1 => Print no title
		if "-vertical" in sys.argv:
			xtitle1 = 0.41
		else:
			xtitle1 = 0.5
		if "-xtitle" in sys.argv:
			xtitle1 = float(sys.argv[sys.argv.index("-xtitle")+1])
		if title2 != '':
			fig.suptitle(title2, y=1.02, x=xtitle1)
		else:
			fig.suptitle(r"Observations", y=1.02, x=xtitle1)

	#########################
	# Set Legend and Limits #
	#########################
	
	plt.savefig(savepath + "stokes_obs" + add)

	##############################################
	#  Plot I, Q, U and V  at wave for result	#
	##############################################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True,
			 gridspec_kw=dict(hspace=0), figsize=figsize)
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if (I2.shape[1] - I2.shape[0]) > 100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0],figsize[1]/4],
										layout="compressed",
									)
		else:
			fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=figsize,
										layout="compressed",
									)
		#fig.subplots_adjust(hspace=0.0)

	############################
	# Plot the Stokes profiles #
	############################
	im1 = ax1.imshow(I2[:,:,wave_ind2]  .transpose(), origin=origin, vmin = limits_stokes2[0][0], vmax = limits_stokes2[0][1], extent=Map_plot)
	im2 = ax2.imshow(Q2[:,:,wave_ind2]  .transpose(), origin=origin, vmin = limits_stokes2[1][0], vmax = limits_stokes2[1][1], cmap = 'PuOr', extent=Map_plot)
	im3 = ax3.imshow(U2[:,:,wave_ind2]  .transpose(), origin=origin, vmin = limits_stokes2[2][0], vmax = limits_stokes2[2][1], cmap = 'PuOr', extent=Map_plot)
	im4 = ax4.imshow(V2[:,:,waveV_ind2].transpose(), origin=origin, vmin = limits_stokes2[3][0], vmax = limits_stokes2[3][1], cmap = 'PuOr', extent=Map_plot)

	#####################
	#	Set labels	#
	#####################
	ax1.set_title(r'$\mathrm{I} / \mathrm{I}_c $ @' + "%.3f" % wave + r" \AA")
	ax2.set_title(r'$\mathrm{Q} / \mathrm{I}_c $')
	ax3.set_title(r'$\mathrm{U} / \mathrm{I}_c $')
	ax4.set_title(r'$\mathrm{V} / \mathrm{I}_c $')
	if waveV != -1:
		ax4.set_title(r'$\mathrm{V} / \mathrm{I}_c$ @' + "%.3f" % waveV + r" $\AA$")		


	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.057 * frac, pad=0.04)
	#cbar1.set_label(label = r'$I / I_c $', loc = 'center')
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.057 * frac, pad=0.04)
	#cbar2.set_label(label = r'$Q / I_c $', loc = 'center')
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.057 * frac, pad=0.04)
	#cbar3.set_label(label = r'$U / I_c $', loc = 'center')
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.057 * frac, pad=0.04)
	#cbar4.set_label(label = r'$V / I_c $', loc = 'center')
	############

	#####################
	#	Set labels	#
	#####################
	if "-vertical" not in sys.argv:
		ax3.set_xlabel(f"x [{units}]")
	ax4.set_xlabel(f"x [{units}]")
	ax1.set_ylabel(f"y [{units}]")
	ax3.set_ylabel(f"y [{units}]")

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################

	if title1 != "-1": # -1 => Print no title
		if "-vertical" in sys.argv:
			xtitle1 = 0.41
		else:
			xtitle1 = 0.5
		if "-xtitle" in sys.argv:
			xtitle1 = float(sys.argv[sys.argv.index("-xtitle")+1])
		if title1 != '':
			fig.suptitle(title1, y=1.02, x=xtitle1)
		else:
			fig.suptitle(r"Inversion Results", y=1.02, x=xtitle1)

	plt.savefig(savepath + "stokes" + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["_____","-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-chi2", "-Bz"]
	labels = ["", r"$T$ [K]", r"$\log P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", r"$B$ [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", r"$z$ [km]", r"$\log P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$",r"$\chi^2$", r"$B \cdot \cos \gamma$ [G]"]
	titles   = ["",r"Temperature", r"Electron Pressure",r"Microturbulence Velocity", r"Magnetic Field",	r"Line-of-Sight Velocity", r"Inclination", r"Azimuth", r"Height", r"Gas Pressure", r"Density$", r"$\chi^2$", r"Line-of-Sight Magnetic Field"]
	cmap = [None,None,None,None,'cividis','seismic','jet','hsv',None,None,None,'gist_gray', None]
	limits = [[None,None],[np.min(models_inv[:,:,1]),np.max(models_inv[:,:,1])],[None,None],[None,None],[None,None],[None,None],[0,180],[0,180],[None, None],[None, None],[None, None],[0,50],[-2000,2000]]
	i = 0

	if "-limitT" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitT")+1].split(',')
		limits[1] = [int(i) for i in temp ]
	
	if "-limitB" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitB")+1].split(',')
		limits[4] = [int(i) for i in temp ]
	if "-limitchi2" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitchi2")+1].split(',')
		limits[11] = [int(i) for i in temp ]
	
	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			# Plot
			fig, ax = plt.subplots(figsize=[figsize[0]/2,figsize[1]/2], layout="compressed")
			ax.set_title(titles[i] + r" @ $\log \tau = $" + str(taus[i]))
			if inputs[i] == "-Bz":
				im = ax.imshow((models_inv.B*np.cos(models_inv.gamma*np.pi/180)).transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			elif inputs[i] == '-chi2':
				im = ax.imshow(chi2_inv.transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			else:
				im = ax.imshow(models_inv.get_attribute(inputs[i][2:]).T, cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)

			# Set labels
			ax.set_xlabel(f"x [{units}]")
			ax.set_ylabel(f"y [{units}]")
			############
			# Colorbar #
			cbar = fig.colorbar(im, ax=ax, fraction=0.057 * frac, pad=0.04)
			cbar.set_label(label = labels[i], loc = 'center')
			############
			# set the spacing between subplots
			#plt.tight_layout(pad=2)
			plt.savefig(savepath + "plot_" + str(inputs[i][1:]) + add)


	# Plot T,B,vlos, inc in one figure
	titles   = ["",r"Temperature", r"Electron Pressure",r"Microturbulence Velocity", r"Magnetic Field",	r"Line-of-Sight Vel.", r"Inclination", r"Azimuth", r"Height", r"Gas Pressure", r"Density$", r"$\chi^2$", r"Line-of-Sight Magnetic Field"]

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True,
			 gridspec_kw=dict(hspace=0), figsize=figsize)
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if (I2.shape[1] - I2.shape[0]) > 100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0],figsize[1]/4],
										layout="compressed",
									)
		else:
			fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=figsize,
										layout="compressed",
									)
		#fig.subplots_adjust(hspace=0.0)


	im1 = ax1.imshow(models_inv.T.transpose(), cmap=cmap[1], origin = origin, vmin = limits[1][0], vmax = limits[1][1],extent=Map_plot)
	im2 = ax2.imshow(models_inv.B.transpose(), cmap=cmap[4], origin = origin, vmin = limits[4][0], vmax = limits[4][1],extent=Map_plot)
	im3 = ax3.imshow(models_inv.vlos.transpose(), cmap=cmap[5], origin = origin, vmin = limits[5][0], vmax = limits[5][1],extent=Map_plot)
	if "-plot_chi2" in sys.argv:
		im4 = ax4.imshow(chi2_inv.transpose(), cmap=cmap[11], origin = origin, vmin = limits[11][0], vmax = limits[11][1],extent=Map_plot)
	else:
		im4 = ax4.imshow(models_inv.gamma.transpose(), cmap=cmap[6], origin = origin, vmin = limits[6][0], vmax = limits[6][1],extent=Map_plot)

	#####################
	#	Set labels	#
	#####################
	#if "-vertical" not in sys.argv:
	#	ax3.set_xlabel(r"x [Pixels]")
	if "-vertical" not in sys.argv:
		ax1.set_xlabel(f"x [{units}]")
		ax2.set_xlabel(f"x [{units}]")
		ax3.set_xlabel(f"x [{units}]")
	ax4.set_xlabel(f"x [{units}]")
	ax1.set_ylabel(f"y [{units}]")
	ax2.set_ylabel(f"y [{units}]")
	ax3.set_ylabel(f"y [{units}]")
	ax4.set_ylabel(f"y [{units}]")

	###############################
	#	Set title and legend	#
	###############################
	if "-vertical" not in sys.argv:
		ax1.set_title(titles[1] + r" @ $\log\tau = $ " + str(logT), fontsize=16)
		ax2.set_title(titles[4] + r" @ $\log\tau = $ " + str(logB), fontsize=16)
		ax3.set_title(titles[5] + r" @ $\log\tau = $ " + str(logV), fontsize=16)
		if "-plot_chi2" in sys.argv:
			ax4.set_title(titles[11], fontsize=16)
		else:
			ax4.set_title(titles[6] + r" @ $\log\tau = $ " + str(logI), fontsize=16)

	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.057 * frac, pad=0.04)
	cbar1.set_label(label = labels[1], loc = 'center')
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.057 * frac, pad=0.04)
	cbar2.set_label(label = labels[4], loc = 'center')
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.057 * frac, pad=0.04)
	cbar3.set_label(label = labels[5], loc = 'center')
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.057 * frac, pad=0.04)
	if "-plot_chi2" in sys.argv:
		cbar4.set_label(label = labels[11], loc = 'center')	
	else:
		cbar4.set_label(label = labels[6], loc = 'center')
	############
	# Set title position depending on the chosen plot and consider the flags hinode and gris
	if title3 != "-1":
		if "-vertical" in sys.argv:
			xtitle1 = 0.41
		else:
			xtitle1 = 0.5
		if "-xtitle" in sys.argv:
			xtitle1 = float(sys.argv[sys.argv.index("-xtitle")+1])
		if title3 != '':
			fig.suptitle(title3, y=1.02, x=xtitle1)
		else:
			fig.suptitle(r"Result Models", y=1.02, x=xtitle1)

	plt.savefig(savepath + "inversion" + add)

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])
	plot(conf, float(sys.argv[2]), float(sys.argv[3]))






