"""
Plots the result of the SIR inversion
"""

import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import os, sys
sys.path.append(sys.path[0] + "/../..")
sys.path.append(sys.path[0] + "/../..")
import sir, obs
import definitions as d
from change_config_path import change_config_path

# Use profile_stk and model

def help():
	"""
	Help Page
	"""
	print("result - Plots the result of the inversion")
	print("Usage: python result.py [OPTION]")
	print()
	print("1:         Config")
	print("2:         wavelength at which stuff is plotted in A")
	print("3:         tau value at which stuff is plotted")
	print()
	print("-data:     Rel. path to the spectral veil corrected data if standard labelling is not used, optional.")	
	print("-stokes:   Rel. path to the Stokes result if standard labelling is not used, optional.")
	print("-models1:  Rel. path to the Models 1 of the inversion if standard labelling is not used.")
	print("-models2:  Rel. path to the Models 2 of the inversion if standard labelling is not used.")
	print("-chi:      Rel. path to the chi2 file of the inversion if standard labelling is not used.")
	print("-save:     Additional save path (optional, default './')")
	print("-add:      Additional text in filenames (optional)")
	print("-label:    Additional text in the title of the single plots")
	print("-title1:   Title in Result Stokes plot")
	print("-title2:   Title in Obs. Stokes plot")
	print("-title3:   Title in Model 1 plot with 4 plots")
	print("-title4:   Title in Model 2 plot with 4 plots")
	print("-T:        Plot temperature in K")
	print("-Pe:       Plot electron pressure in dyn/cm^2")
	print("-vmicro:   Plot microturbulence in cm/s")
	print("-B:        Plot magentic field strength in Gauss")
	print("-vlos:     Plot line of sight velocity in km/s")
	print("-inc:      Plot inclination by subtracting in deg")
	print("-azi:      Plot azimuth by adding in deg")
	print("-z:        Plot real height in km")
	print("-Pg:       Plot gas pressure in dyn/cm^2")
	print("-rho:      Plot density")
	print("-vertical: Plot spectra vertically")
	print("-chi2:     Plot chi2")
	print("-fill:     Plot the filling factor")
	print("-plot_chi2: Plot chi2 in the 4 subplots plot")
	print("-plot_fill: Plot filling factor in the 4 subplots plot")
	print("-waveQ:    Plot Stokes Q in another wavelength position in A")
	print("-waveU:    Plot Stokes U in another wavelength position in A")
	print("-waveV:    Plot Stokes V in another wavelength position in A")
	print("-logT:     Plot Temperature at this log tau.")
	print("-logB:     Plot Magnetic Field at this log tau.")
	print("-logV:     Plot LoS Velocity at this log tau.")
	print("-logI:     Plot Inclination at this log tau.")
	print("-limitxy:  Limit in the x and y plot as a list xmin,xmax,ymin,xmax")
	print("-limitT:   Set the limit for the colorbar in the temperature.")
	print("-limitB:   Set the limit for the colorbar in the magnetic field.")
	print("-limitchi2:Set the limit for the colorbar in chi2.")
	print("-limitI:   Set the limit for the colorbar in Stokes I.")
	print("-arc:      Print x and y axis in arcseconds")
	print("-flipx:    Mirror/Flip data as sometimes it is wrong in GRIS with the location on the sun")

	sys.exit()

def check_range(range_wave_ang, wave):
	"""
	Check if the given wavelength is in the range

	Parameter
	---------
	range_wave_ang : numpy nd array
		Array containing the wavelength ranges in A for the given inversion. The size is nx2 depending on the Grid file
	wave : float
		Wavelength in A to be checked

	Return
	------
	Wavelength in the range

	"""
	# Check for the wavelength if it is in the inversion range:
	if wave < np.min(range_wave_ang):
		print("Note that the given wavelength is not in the inversion range. Take closest one in range...")
		wave = np.min(range_wave_ang)
		return wave
	elif wave > np.max(range_wave_ang):
		print("Note that the given wavelength is not in the inversion range. Take closest one in range...")
		wave = np.max(range_wave_ang)
		return wave

	in_range = False
	for i in range_wave_ang:
		# Check that the wavelength is in the range
		if i[0] <= wave <= i[1]:
			return wave
	if not in_range:
		print("Note that the given wavelength is not in the inversion range. Take closest one in range...")
		wave = wave.flatten()[np.argmin(abs(range_wave_ang-wave))]
	return wave



def plot(conf, wave, tau, Type = "_1", plot_stokes = True):
	"""
	Plots the result of the inversion

	Parameter
	---------
	conf : dict
		Dict. with all the information from the config
	wave : float
		Wavelength in A where the Stokes vector is plottet
	tau : float
		log tau value where the model is plotted
	Type : string, optional
		prefix for determining which model is used. Default: '_1'
	plot_stokes : bool, optional
		Plot the stokes vector. Default: True

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
			change_config_path(conf,os.path.abspath(os.getcwd()))

	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	path = conf["path"]
	waves = np.load(os.path.join(path, conf['waves']))
	Map = conf['map']
	range_wave_ang = sir.pixel_to_angstrom(waves, conf['range_wave'])
	
	# Check whether the wavelength is in range
	wave = check_range(range_wave_ang, wave)
	
	if "-waveQ" in sys.argv:
		waveQ = float(sys.argv[sys.argv.index("-waveQ")+1])
		waveQ = check_range(range_wave_ang, waveQ)
	else:
		waveQ = wave

	if "-waveU" in sys.argv:
		waveU = float(sys.argv[sys.argv.index("-waveU")+1])
		waveU = check_range(range_wave_ang, waveU)
	else:
		waveU = wave

	if "-waveV" in sys.argv:
		waveV = float(sys.argv[sys.argv.index("-waveV")+1])
		waveV = check_range(range_wave_ang, waveV)
	else:
		waveV = wave
	
	if "-data" not in sys.argv:
		stokes = obs.load_data(conf, filename=conf['cube_inv'])
	else:
		filename = sys.argv[sys.argv.index("-data")+1]
		stokes = obs.load_data(conf, filename=filename)

	if "-stokes" not in sys.argv:
		stokes_inv = np.load(os.path.join(path, conf['inv_out'] + d.end_stokes))
	else:
		filename = sys.argv[sys.argv.index("-stokes")+1]
		stokes_inv = np.load(filename)

	if Type == "_1":
		if "-models1" not in sys.argv:
			models_inv = np.load(os.path.join(path, conf['inv_out'] + d.end_models1))
		else:
			filename = sys.argv[sys.argv.index("-models1")+1]
			models_inv = np.load(filename)
	if Type == "_2":
		if "-models2" not in sys.argv:
			models_inv = np.load(os.path.join(path, conf['inv_out'] + d.end_models2))
		else:
			filename = sys.argv[sys.argv.index("-models2")+1]
			models_inv = np.load(filename)

	if Type == "_1":
		if "-models1" not in sys.argv:
			errors_inv = np.load(os.path.join(path,conf['inv_out'] + d.end_errors1))
		else:
			filename = sys.argv[sys.argv.index("-models1")+1].replace(".mod",".err")
			errors_inv = np.load(filename)
	if Type == "_2":
		if "-models2" not in sys.argv:
			errors_inv = np.load(os.path.join(path,conf['inv_out'] + d.end_errors2))
		else:
			filename = sys.argv[sys.argv.index("-models2")+1].replace(".mod",".err")
			errors_inv = np.load(filename)

	if "-chi" not in sys.argv:
		chi2_inv = np.load(os.path.join(path,conf['chi2']))
	else:
		filename = sys.argv[sys.argv.index("-chi")+1]
		chi2_inv = np.load(filename)

	if "-fill" in sys.argv:
		if Type == "_1":
			fill_inv = np.load(os.path.join(path,d.filling_factor + d.end_models1))
		if Type == "_2":
			fill_inv = np.load(os.path.join(path,d.filling_factor + d.end_models2))

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

	# Cut data in x and y position	
	if "-limitxy" in sys.argv:
		limit_xy = np.array([int(i) for i in sys.argv[sys.argv.index("-limitxy")+1].split(",")], dtype=int)

		# Cut data to the new range:
		stokes_inv = stokes_inv[limit_xy[0]-Map[0]:limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2]:limit_xy[3]+1-(Map[3]+1)]
		models_inv = models_inv[limit_xy[0]-Map[0]:limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2]:limit_xy[3]+1-(Map[3]+1)]
		errors_inv = errors_inv[limit_xy[0]-Map[0]:limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2]:limit_xy[3]+1-(Map[3]+1)]

		# Save the new limits as the new Map
		Map = limit_xy

	instrument = conf["instrument"]  # Instrument used for labels

	# Get waves indeces used in the observation closest to the chosen one
	wave_ind = np.argmin(abs(waves-wave))
	wave = waves[wave_ind]

	waveQ_ind = np.argmin(abs(waves-waveQ))
	waveQ = waves[waveQ_ind]

	waveU_ind = np.argmin(abs(waves-waveU))
	waveU = waves[waveU_ind]

	waveV_ind = np.argmin(abs(waves-waveV))
	waveV = waves[waveV_ind]

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = path + "/" + sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath):
			os.mkdir(savepath)

	# Additional text in output
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]
	
	# Additional text in the title of the single plots
	add_label = ''
	if '-label' in sys.argv:
		add_label = sys.argv[sys.argv.index("-label")+1]

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
	title4 = ''
	if '-title4' in sys.argv:
		title4 = sys.argv[sys.argv.index("-title4")+1]

	# Restrict models to the given tau
	ind = np.argmin(abs(models_inv[0, 0, 0, :] - tau))
	indT = np.argmin(abs(models_inv[0, 0, 0, :] - logT))
	indB = np.argmin(abs(models_inv[0, 0, 0, :] - logB))
	indV = np.argmin(abs(models_inv[0, 0, 0, :] - logV))
	indI = np.argmin(abs(models_inv[0, 0, 0, :] - logI))

	# Cut models to the specific log tau value
	models_inv_temp = np.copy(models_inv)
	models_inv = np.zeros(shape=(models_inv.shape[0],models_inv.shape[1],models_inv.shape[2]))
	for i in [0, 2, 3, 7, 8, 9, 10]:
		models_inv[:, :, i] = models_inv_temp[:, :, i, ind]

	models_inv[:, :, 1] = models_inv_temp[:, :, 1, indT]
	models_inv[:, :, 4] = models_inv_temp[:, :, 4, indB]
	models_inv[:, :, 5] = models_inv_temp[:, :, 5, indV]
	models_inv[:, :, 6] = models_inv_temp[:, :, 6, indI]

	taus = [tau for i in range(11)]
	taus[1] = logT
	taus[4] = logB
	taus[5] = logV
	taus[6] = logI
	
	# Pressure in log
	models_inv[:, :, 2] = np.log(models_inv[:, :, 2])
	models_inv[:, :, 9] = np.log(models_inv[:, :, 9])

	# vlos in km/s
	models_inv[:, :, 5] = models_inv[:, :, 5] / 1e5

	#########################
	#  PLOTTING STUFF		#
	#########################
	if I2.shape[0] / I2.shape[1] >= 0.8:
		alpha = I2.shape[0] / 12
		figsize = [I2.shape[0] / alpha , I2.shape[1]/alpha-1]
		frac = figsize[1] / figsize[0]
	else:
		alpha = I2.shape[1] / 18
		figsize = [I2.shape[0] / alpha, I2.shape[1]/alpha-5]
		frac = figsize[1] / figsize[0] * 1.25

	####################################
	#	Plot settings for arcsecond	#
	####################################
	if "-arc" in sys.argv:
		infos = dict(np.loadtxt(os.path.join(path, d.header_infos), delimiter='=', dtype=str))
		if conf['instrument'] == 'GRIS':
			# y position starts at the end but I start with the pixels lower left
			y_pos = float(infos['CRVAL2']) + (stokes.shape[1] - 1) * float(infos['CDELT2'])
			y_max = float(infos['CRVAL2'])  # Used if flipx is used
			Map_plot = [float(infos['CRVAL1']) - float(infos['CDELT1']) * (Map[0] - 1),
						float(infos['CRVAL1']) - float(infos['CDELT1']) * (Map[1] - 1),
						y_pos - float(infos['CDELT2']) * (Map[2] - 1),
						y_pos - float(infos['CDELT2']) * (Map[3] - 1)
						]
		elif conf['instrument'] == 'Hinode':
			delta_y = float(infos['CDELT2'])  # Delta y of slit
			delta_x = float(infos['XSCALE'])
			x_pos = float(infos['XCEN'])
			y_pos = float(infos['YCEN'])
			y_max = y_pos + ((stokes.shape[1] - 1) - infos['CRPIX2']) * delta_y  # Used for flipx
			Map_plot = [
				(Map[0] - 1) * delta_x + x_pos,
				(Map[1] - 1) * delta_x + x_pos,
				(Map[2] - 1 - infos['CRPIX2']) * delta_y + y_pos,
				(Map[3] - 1 - infos['CRPIX2']) * delta_y + y_pos
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

	if plot_stokes:
		I1 = stokes[Map[0]:Map[1]+1, Map[2]:Map[3]+1, 0, :]
		Q1 = stokes[Map[0]:Map[1]+1, Map[2]:Map[3]+1, 1, :]
		U1 = stokes[Map[0]:Map[1]+1, Map[2]:Map[3]+1, 2, :]
		V1 = stokes[Map[0]:Map[1]+1, Map[2]:Map[3]+1, 3, :]

		I2 = stokes_inv[:, :,0 ,:]
		Q2 = stokes_inv[:, :,1 ,:]
		U2 = stokes_inv[:, :,2 ,:]
		V2 = stokes_inv[:, :,3 ,:]

		limits_stokes1 = [  [np.min(I1[:,:,wave_ind]), np.max(I1[:,:,wave_ind])],
						[-np.max(np.abs(Q1[:,:,wave_ind])), np.max(np.abs(Q1[:,:,waveQ_ind]))],
						[-np.max(np.abs(U1[:,:,wave_ind])), np.max(np.abs(U1[:,:,waveU_ind]))],
						[-np.max(np.abs(V1[:,:,wave_ind])), np.max(np.abs(V1[:,:,waveV_ind]))]
			]
		limits_stokes2 = [  [np.min(I2[:,:,wave_ind]), np.max(I2[:,:,wave_ind])],
						[-np.max(np.abs(Q2[:,:,wave_ind])), np.max(np.abs(Q2[:,:,waveQ_ind]))],
						[-np.max(np.abs(U2[:,:,wave_ind])), np.max(np.abs(U2[:,:,waveU_ind]))],
						[-np.max(np.abs(V2[:,:,wave_ind])), np.max(np.abs(V2[:,:,waveV_ind]))]
			]


		if "-limitI" in sys.argv:
			temp = sys.argv[sys.argv.index("-limitI")+1].split(',')
			limits_stokes1[0] = [int(i) for i in temp ]
			limits_stokes2[0] = [int(i) for i in temp ]

		##############################################
		#  Plot I, Q, U and V  at wave for obs		#
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



		############################
		# Plot the Stokes profiles #
		############################
		im1 = ax1.imshow(I1[:,:,wave_ind]  .transpose(), origin=origin, vmin = limits_stokes1[0][0], vmax = limits_stokes1[0][1], extent=Map_plot)
		im2 = ax2.imshow(Q1[:,:,waveQ_ind]  .transpose(), origin=origin, vmin = limits_stokes1[1][0], vmax = limits_stokes1[1][1], cmap = 'PuOr', extent=Map_plot)
		im3 = ax3.imshow(U1[:,:,waveU_ind]  .transpose(), origin=origin, vmin = limits_stokes1[2][0], vmax = limits_stokes1[2][1], cmap = 'PuOr', extent=Map_plot)
		im4 = ax4.imshow(V1[:,:,waveV_ind].transpose(), origin=origin, vmin = limits_stokes1[3][0], vmax = limits_stokes1[3][1], cmap = 'PuOr', extent=Map_plot)

		#####################
		#	Set labels	#
		#####################
		ax1.set_title(r'$I / I_c $ @' + "%.3f" % wave + r" \AA")
		ax2.set_title(r'$Q / I_c $')
		ax3.set_title(r'$U / I_c $')
		ax4.set_title(r'$V / I_c $')

		
		if waveQ != wave:
			ax2.set_title(r'$Q / I_c$ @' + "%.3f" % waveQ + r" \AA")	
		if waveU != wave:
			ax3.set_title(r'$U / I_c$ @' + "%.3f" % waveU + r" \AA")	
		if waveV != wave:
			ax4.set_title(r'$V / I_c$ @' + "%.3f" % waveV + r" \AA")	


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
			if title2 != '':
				fig.suptitle(title2, y=1.02, x=xtitle1)
			else:
				fig.suptitle(r"Observation Stokes Vector", y=1.02, x=xtitle1)


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


		############################
		# Plot the Stokes profiles #
		############################
		im1 = ax1.imshow(I2[:,:,wave_ind ]  .transpose(), origin=origin, vmin = limits_stokes2[0][0], vmax = limits_stokes2[0][1], extent=Map_plot)
		im2 = ax2.imshow(Q2[:,:,waveQ_ind]  .transpose(), origin=origin, vmin = limits_stokes2[1][0], vmax = limits_stokes2[1][1], cmap = 'PuOr', extent=Map_plot)
		im3 = ax3.imshow(U2[:,:,waveU_ind]  .transpose(), origin=origin, vmin = limits_stokes2[2][0], vmax = limits_stokes2[2][1], cmap = 'PuOr', extent=Map_plot)
		im4 = ax4.imshow(V2[:,:,waveV_ind]  .transpose(), origin=origin, vmin = limits_stokes2[3][0], vmax = limits_stokes2[3][1], cmap = 'PuOr', extent=Map_plot)

		#####################
		#	Set labels	#
		#####################
		ax1.set_title(r'$I / I_c $ @' + "%.3f" % wave + r" \AA")
		ax2.set_title(r'$Q / I_c $')
		ax3.set_title(r'$U / I_c $')
		ax4.set_title(r'$V / I_c $')
		if waveQ != wave:
			ax2.set_title(r'$Q / I_c$ @' + "%.3f" % waveQ + r" \AA")
		if waveU != wave:
			ax3.set_title(r'$U / I_c$ @' + "%.3f" % waveU + r" \AA")
		if waveV != wave:
			ax4.set_title(r'$V / I_c$ @' + "%.3f" % waveV + r" \AA")


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
			if title1 != '':
				fig.suptitle(title1, y=1.02, x=xtitle1)
			else:
				fig.suptitle(r"Inversion Results", y=1.02, x=xtitle1)

		plt.savefig(savepath + "stokes" + add)

	###################################################
	#			Plot physical parameters			#
	###################################################

	# Define labels and get from arguments which parameter should be plot
	inputs = ["_____","-T", '-Pe', '-vmicro', '-B', "-vlos", "-inc", "-azi", "-z", "-Pg","-rho","-chi2", "-Bz", "-fill"]
	labels = ["", r"$T$ [K]", r"$\log P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", r"$B$ [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", r"$z$ [km]", r"$\log P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$",r"$\chi^2$", r"$B \cdot \cos \gamma$ [G]", r"$\alpha$"]
	titles   = ["",r"Temperature", r"Electron Pressure",r"Microturbulence Velocity", r"Magnetic Field",	r"Line-of-Sight Velocity", r"Inclination", r"Azimuth", r"Height", r"Gas Pressure", r"Density$", r"$\chi^2$", r"Line-of-Sight Magnetic Field", r"Filling Factor"]
	cmap = [None,None,None,None,'cividis','seismic','jet','hsv',None,None,None,'gist_gray', None, None]
	limits = [[None,None],[np.min(models_inv[:,:,1]),np.max(models_inv[:,:,1])],[None,None],[None,None],[0,4000],[-5,5],[0,180],[0,180],[None, None],[None, None],[None, None],[0,50],[-2000,2000], [0,1]]
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
				
			if inputs[i] == "-chi2":
				ax.set_title(titles[i]  + add_label)
			elif inputs[i] == "-fill":
				ax.set_title(titles[i] + f" for Model {Type[1]}" + add_label)			

			else:
				ax.set_title(titles[i] + r" @ $\log \tau = $" + str(taus[i]) + f" for Model {Type[1]}" + add_label)
			
			if inputs[i] == "-Bz":
				im = ax.imshow((models_inv[:,:,4]*np.cos(models_inv[:,:,6]*np.pi/180)).transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			elif inputs[i] == '-chi2':
				im = ax.imshow(chi2_inv.transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			elif inputs[i] == '-fill':
				im = ax.imshow(fill_inv.transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			else:
				im = ax.imshow(models_inv[:,:,i], cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)

			# Set labels
			ax.set_xlabel(f"x [{units}]")
			ax.set_ylabel(f"y [{units}]")
			############
			# Colorbar #
			cbar = fig.colorbar(im, ax=ax, fraction=0.057 * frac, pad=0.04)
			cbar.set_label(label = labels[i], loc = 'center')

			plt.savefig(savepath + "plot_" + str(inputs[i][1:]) + Type + add)


	# Plot T,B,vlos, inc in one figure
	titles   = ["",r"Temperature", r"Electron Pressure",r"Microturbulence Velocity", r"Magnetic Field",	r"Line-of-Sight Vel.", r"Inclination", r"Azimuth", r"Height", r"Gas Pressure", r"Density$", r"$\chi^2$", r"Line-of-Sight Magnetic Field", r"Filling Factor"]

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


	im1 = ax1.imshow(models_inv[:,:,1].transpose(), cmap=cmap[1], origin = origin, vmin = limits[1][0], vmax = limits[1][1],extent=Map_plot)
	im2 = ax2.imshow(models_inv[:,:,4].transpose(), cmap=cmap[4], origin = origin, vmin = limits[4][0], vmax = limits[4][1],extent=Map_plot)
	im3 = ax3.imshow(models_inv[:,:,5].transpose(), cmap=cmap[5], origin = origin, vmin = limits[5][0], vmax = limits[5][1],extent=Map_plot)
	if "-plot_chi2" in sys.argv:
		im4 = ax4.imshow(chi2_inv.transpose(), cmap=cmap[11], origin = origin, vmin = limits[11][0], vmax = limits[11][1],extent=Map_plot)
	elif "-plot_fill" in sys.argv:
		im4 = ax4.imshow(fill_inv.transpose(), cmap=cmap[13], origin = origin, vmin = limits[13][0], vmax = limits[13][1],extent=Map_plot)
	else:
		im4 = ax4.imshow(models_inv[:,:,6].transpose(), cmap=cmap[6], origin = origin, vmin = limits[6][0], vmax = limits[6][1],extent=Map_plot)

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
		elif "-plot_fill" in sys.argv:
			ax4.set_title(titles[13], fontsize=16)
		else:
			ax4.set_title(titles[6] + r" @ $\log\tau = $ " + str(logI), fontsize=16)

	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.046 * frac, pad=0.04)
	cbar1.set_label(label = labels[1], loc = 'center')
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046 * frac, pad=0.04)
	cbar2.set_label(label = labels[4], loc = 'center')
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.046 * frac, pad=0.04)
	cbar3.set_label(label = labels[5], loc = 'center')
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.046 * frac, pad=0.04)
	if "-plot_chi2" in sys.argv:
		cbar4.set_label(label = labels[11], loc = 'center')
	elif "-plot_fill" in sys.argv:
		cbar4.set_label(label = labels[13], loc = 'center')		
	else:
		cbar4.set_label(label = labels[6], loc = 'center')
	############
	# Set title position depending on the chosen plot and consider the flags hinode and gris
	if Type == "_1":
		if title3 != "-1":
			if "-vertical" in sys.argv:
				xtitle1 = 0.41
			else:
				xtitle1 = 0.55
			if title3 != '':
				fig.suptitle(title3, y=1.02, x=xtitle1)
			else:
				fig.suptitle(f"Inversion Results for Model {Type[1]}", y=1.02, x=xtitle1)
	else:
		if title4 != "-1":
			if "-vertical" in sys.argv:
				xtitle1 = 0.41
			else:
				xtitle1 = 0.5
			if title4 != '':
				fig.suptitle(title4, y=1.02, x=xtitle1)
			else:
				fig.suptitle(f"Inversion Results for Model {Type[1]}", y=1.02, x=xtitle1)

	plt.savefig(savepath + "inversion" + Type + add)

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])
	plot(conf, float(sys.argv[2]), float(sys.argv[3]))
	plot(conf, float(sys.argv[2]), float(sys.argv[3]), Type = "_2" , plot_stokes = False)






