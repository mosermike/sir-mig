"""

Plots the result of the SIR inversion

"""

import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
import os
import sys
sys.path.append(sys.path[0] + "/../src")
import sir
import definitions as d
import model_atm as m
import profile_stk as p
import chi2_stk as c

# TODO do the figure size as in gris_firtez

def _help():
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
	sir.option("-models: (1C)","Rel. path to the Models of the inversion if standard labelling is not used.")
	sir.option("-models1: (2C)","Rel. path to the Models 1 of the inversion if standard labelling is not used.")
	sir.option("-models2: (2C)","Rel. path to the Models 2 of the inversion if standard labelling is not used.")
	#sir.option("-errors:","Rel. path to the Errors of the inversion if standard labelling is not used.")
	sir.option("-chi:","Rel. path to the chi2 file of the inversion if standard labelling is not used.")
	sir.option("-save:","Additional save path (optional, default './')")
	sir.option("-add:","Additional text in filenames (optional)")
	sir.option("-label:","Add label text (optional)")
	sir.option("-title1:","Title in Result Stokes plot")
	sir.option("-title2:","Title in Obs. Stokes plot")
	sir.option("-title3: (1C, 2C)","Title in Model 1 plot")
	sir.option("-title4: (2C)","Title in Model 2 plot with 4 plots")
	sir.option("-xtitle:","Changing the x position of the title")
	sir.option("-T:","Plot temperature in K")
	sir.option("-Pe:","Plot electron pressure in dyn/cm^2")
	sir.option("-vmicro:","Plot microturbulence in cm/s")
	sir.option("-B:","Plot magentic field strength in Gauss")
	sir.option("-vlos:","Plot line of sight velocity in km/s")
	sir.option("-gamma:","Plot inclination by subtracting in deg")
	sir.option("-phi:","Plot azimuth by adding in deg")
	sir.option("-z:","Plot real height in km")
	sir.option("-Pg:","Plot gas pressure in dyn/cm^2")
	sir.option("-rho:","Plot density")
	sir.option("-vertical:","Plot spectra vertically")
	sir.option("-chi2:","Plot chi2")
	sir.option("-fill: (2C)","Plot the filling factor")
	sir.option("-plot_chi2:","Plot chi2 in the 4 subplots plot")
	sir.option("-plot_fill: (2C)","Plot filling factor in the 4 subplots plot")
	sir.option("-waveQ:","Plot Stokes Q in another wavelength position in A")
	sir.option("-waveU:","Plot Stokes U in another wavelength position in A")
	sir.option("-waveV:","Plot Stokes V in another wavelength position in A")
	sir.option("-logT:","Plot Temperature at this log tau.")
	sir.option("-logB:","Plot Magnetic Field at this log tau.")
	sir.option("-logV:","Plot LoS Velocity at this log tau.")
	sir.option("-logI:","Plot Inclination at this log tau.")
	sir.option("-limitxy:","Limit in the x and y plot as a list xmin,xmax,ymin,xmax")
	sir.option("-limitT:","Set the limit for the colorbar in the temperature.")
	sir.option("-limitB:","Set the limit for the colorbar in the magnetic field.")
	sir.option("-limitv:","Set the limit for the colorbar in the velocity.")
	sir.option("-limitchi2:","Set the limit for the colorbar in chi2.")
	sir.option("-limitI:","Set the limit for the colorbar in Stokes I.")
	sir.option("-arc:","Print x and y axis in arcseconds")
	sir.option("-kG:","Plot in kG (Kilogauss)")
	sir.option("-flipy:","Flip image along y axis as sometimes the orientation is wrong in GRIS with the location on the sun")
	sir.option("-f [float]","Factor for the fontsizes")
	sir.option("-symv","Symmetric limits for vlos")
	sir.option("-rot90","Rotate the image 90 deg")
	sys.exit()


def _check_range(wave_inv, wave):
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
		return wave
	elif wave > np.max(wave_inv):
		print("Note that the given wavelength is not in the inversion range. Take closest one in range...")
		wave = np.max(wave_inv)
		return wave

	return wave


def plot_chi2(figsize, frac, chi2, Map_plot, units, savepath, add, origin, title4 =""):
	"""
	Plots the models

	Parameters
	----------
	figsize : Array
		Array with figure size
	frac : float
		Fraction between y and x
	models_inv : class
		Instance of class model from firtez_dz
	zs : array
		Array containing the indexes in z
	Map_plot : array
		Ticks of x and y
	units : string
		Unit in arcsec or Pixel
	savepath : string
		savepath for the files
	add : string
		Additional text in the saved files
	origin : string
		Orientation of the plots
	title4 : string (not implemented)
		Title for the 4 plots
	
	"""
	
	#######################
	#  Plot I, Q, U and V #
	#######################
	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True,
			 gridspec_kw=dict(hspace=0), figsize=figsize)
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if (chi2.stki.shape[1] - chi2.stki.shape[0]) >= -100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0]*2,figsize[1]*2/4],
										layout="compressed",
									)
		else:
			fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=figsize,
										layout="compressed",
									)
		#fig.subplots_adjust(hspace=0.0)
		
	vmaxs = [None, None, None, None, None]
	cmap = 'Greys'
	ext_cbar = 'neither'
	if "-limitchi2" in sys.argv:
		vmaxs = [ float(i) for i in sys.argv[sys.argv.index("-limitchi2")+1].split(',')]
		import matplotlib as mpl
		cmap = mpl.colormaps[cmap]
		cmap.set_extremes(over='yellow')
		ext_cbar='max'

	############################
	# Plot the Stokes profiles #
	############################
	im1 = ax1.imshow(chi2.stki.transpose(), origin=origin, cmap = cmap, extent=Map_plot, vmax = vmaxs[1])
	im2 = ax2.imshow(chi2.stkq.transpose(), origin=origin, cmap = cmap, extent=Map_plot, vmax = vmaxs[2])
	im3 = ax3.imshow(chi2.stku.transpose(), origin=origin, cmap = cmap, extent=Map_plot, vmax = vmaxs[3])
	im4 = ax4.imshow(chi2.stkv.transpose(), origin=origin, cmap = cmap, extent=Map_plot, vmax = vmaxs[4])

	#####################
	#	Set labels	#
	#####################
	ax1.set_title(r"$\chi^2$ for Stokes $I$")
	ax2.set_title(r"$\chi^2$ for Stokes $Q$")
	ax3.set_title(r"$\chi^2$ for Stokes $U$")
	ax4.set_title(r"$\chi^2$ for Stokes $V$")


	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.057 * frac, pad=0.04, cmap=cmap, extend=ext_cbar)
	#cbar1.set_label(label = r'$I / I_c $', loc = 'center')
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.057 * frac, pad=0.04, cmap=cmap, extend=ext_cbar)
	#cbar2.set_label(label = r'$Q / I_c $', loc = 'center')
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.057 * frac, pad=0.04, cmap=cmap, extend=ext_cbar)
	#cbar3.set_label(label = r'$U / I_c $', loc = 'center')
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.057 * frac, pad=0.04, cmap=cmap, extend=ext_cbar)
	#cbar4.set_label(label = r'$V / I_c $', loc = 'center')
	############

	#####################
	#	Set labels	#
	#####################
	if "-vertical" not in sys.argv:
		ax1.set_xlabel(f"x [{units}]")
		ax2.set_xlabel(f"x [{units}]")
		ax3.set_xlabel(f"x [{units}]")
		ax2.set_ylabel(f"y [{units}]")
		ax4.set_ylabel(f"y [{units}]")
	ax4.set_xlabel(f"x [{units}]")
	ax1.set_ylabel(f"y [{units}]")
	ax3.set_ylabel(f"y [{units}]")

	##################################################################
	# Set title											#
	# The position is relative to the chosen plot (vertical or not)  #
	##################################################################
	if title4 != "-1": # -1 => Print no title
		if "-vertical" in sys.argv:
			xtitle1 = 0.41
		else:
			xtitle1 = 0.5
		if title4 != '':
			fig.suptitle(title4, y=1.02, x=xtitle1)


	plt.savefig(savepath + "chi2_stokes" + add)

	# Plot
	fig, ax = plt.subplots(figsize=[figsize[0]/2,figsize[1]/2], layout="compressed")
	#fig.subplots_adjust(hspace=0.0)
	ax.set_title(r"$\chi^2_{tot}$")
	im = ax.imshow(chi2.tot.transpose(), cmap=cmap, origin = origin,vmax = vmaxs[0],# vmin = limits[i][0], vmax = limits[i][1],
						extent=Map_plot)
	# Set labels
	ax.set_xlabel(f"x [{units}]")
	ax.set_ylabel(f"y [{units}]")
	############
	# Colorbar #
	cbar = fig.colorbar(im, ax=ax, fraction=0.057 * frac, pad=0.04, cmap=cmap, extend=ext_cbar)
	cbar.set_label(label = r"$\chi^2$", loc = 'center')
	############
	plt.savefig(savepath + "chi2_total" + add)
	

def _plot_model(models_inv, tau, figsize, frac, units, title3, title4, savepath, add, chi2, Map_plot, origin, sign1, sign2, n, dx, dy, Type=""):


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

	taus = [tau for i in range(11)]
	taus[1] = logT
	taus[4] = logB
	taus[5] = logV
	taus[6] = logI

	# Restrict models to the given tau
	ind  = np.argmin(abs(models_inv.tau - tau))
	indT = np.argmin(abs(models_inv.tau - logT))
	indB = np.argmin(abs(models_inv.tau - logB))
	indV = np.argmin(abs(models_inv.tau - logV))
	indI = np.argmin(abs(models_inv.tau - logI))

	# Correct for 180 deg ambiguity
	models_inv.correct_phi()

	# Cut models to the specific log tau value
	models_inv.nval = 1
	models_inv.T = models_inv.T[:,:,indT]
	models_inv.Pe = models_inv.Pe[:,:,ind]
	models_inv.vmicro = models_inv.vmicro[:,:,ind]
	models_inv.B = models_inv.B[:,:,indB]
	models_inv.vlos = models_inv.vlos[:,:,indV]
	models_inv.gamma = models_inv.gamma[:,:,indI]
	models_inv.phi = models_inv.phi[:,:,ind]
	models_inv.z = models_inv.z[:,:,ind]
	models_inv.rho = models_inv.rho[:,:,ind]
	models_inv.Pg = models_inv.Pg[:,:,ind]

	# Pressure in log
	models_inv.Pe = np.log(models_inv.Pe)
	models_inv.Pg = np.log(models_inv.Pg)

	# in kG
	if "-kG" in sys.argv:
		models_inv.B /= 1e3


	# Define labels and get from arguments which parameter should be plot
	inputs = ["_____","-T", '-Pe', '-vmicro', '-B', "-vlos", "-gamma", "-phi", "-z", "-Pg","-rho","-Bz","-fill"]
	labels = ["", r"$T$ [K]", r"$\log P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", r"$B$ [G]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", r"$z$ [km]", r"$\log P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$", r"$B$ [G]",r"$\alpha$"]
	titles   = ["",r"Temperature", r"Electron Pressure",r"Microturbulence Velocity", r"Magnetic Field",	r"Line-of-Sight Velocity", r"Inclination", r"Azimuth", r"Height", r"Gas Pressure", r"Density$", r"Magnetic Field $B_{\text{los}}$","Filling Factor"]
	cmap = [None,None,None,None,'cividis','seismic','jet','hsv',None,None,None,None, "gist_gray"]
	limits = [[None,None],[np.min(models_inv.T),np.max(models_inv.T)],[None,None],[None,None],
		   [None,None],[None,None],[0,180],[0,180],[None, None],[None, None],[None, None],[-np.max(np.abs(models_inv.B*np.cos(models_inv.gamma/180*np.pi))),np.max(np.abs(models_inv.B*np.cos(models_inv.gamma/180*np.pi)))], [0,1]]
	i = 0
	if "-kG" in sys.argv:
		labels = ["", r"$T$ [K]", r"$\log P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\mathrm{v}_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", r"$B$ [kG]", r"$\mathrm{v}_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]", r"$z$ [km]", r"$\log P_g$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$\rho$ $\left[\mathrm{dyn}\mathrm{cm}^{-3}\right]$", r"$B$ [kG]",r"$\alpha$"]
	if "-symv" in sys.argv:
		limits[5] = [-np.max(np.abs(models_inv.vlos)),np.max(np.abs(models_inv.vlos))]

	if "-limitT" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitT")+1].split(',')
		limits[1] = [int(i) for i in temp ]
	
	if "-limitB" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitB")+1].split(',')
		limits[4] = [int(i) for i in temp ]
	if "-limitv" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitv")+1].split(',')
		limits[5] = [int(i) for i in temp ]
	if "-limitchi2" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitchi2")+1].split(',')
		limits[11] = [int(i) for i in temp ]
	if Type != "":
		add_str =  f" for Model {Type} "
	else:
		add_str = ""
	for i in range(len(inputs)):
		if inputs[i] in sys.argv:
			# Plot
			fig, ax = plt.subplots(figsize=[figsize[0]/2,figsize[1]/2], layout="compressed")
			if inputs[i] == "-Bz":
				ax.set_title(titles[i] + add_str + r" @ $\log \tau = $" + str(taus[4]))
			elif inputs[i] == "-fill":
				ax.set_title(titles[i] + add_str)
			else:
				ax.set_title(titles[i] + add_str + r" @ $\log \tau = $" + str(taus[i]))
			
			if inputs[i] == "-Bz":
				im = ax.imshow((models_inv.B*np.cos(models_inv.gamma*np.pi/180)).transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			elif inputs[i] == "-fill":
				im = ax.imshow(models_inv.fill.transpose(), cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			else:
				im = ax.imshow(models_inv.get_attribute(inputs[i][1:]).T, cmap=cmap[i], origin = origin, vmin = limits[i][0], vmax = limits[i][1],
							extent=Map_plot)
			if inputs[i] == "-vlos":
				if "-arc" in sys.argv and "-varrow" in sys.argv:
					l = np.min([Map_plot[1],Map_plot[3]])*0.025
					ax.arrow(Map_plot[1]/2, Map_plot[3]/2, sign1*abs(dx)*n, sign2*abs(dy)*n, head_width=l, head_length=l, fc='black', ec='black')
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
			plt.savefig(savepath + "plot_" + str(inputs[i][1:]) + Type + add)


	# Plot T,B,vlos, inc in one figure
	titles   = ["",r"Temperature", r"Electron Pressure",r"Microturbulence Velocity", r"Magnetic Field",	r"Line-of-Sight Vel.", r"Inclination", r"Azimuth", r"Height", r"Gas Pressure", r"Density$", r"$\chi^2$", r"Magnetic Field $B_{\text{los}}$"]

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True,
			 gridspec_kw=dict(hspace=0), figsize=figsize)
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if (models_inv.ny - models_inv.nx) >= -100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0]*2,figsize[1]*2/4],
										layout="compressed",
									)
		else:
			fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=figsize,
										layout="compressed",
									)



	im1 = ax1.imshow(models_inv.T.transpose(), cmap=cmap[1], origin = origin, vmin = limits[1][0], vmax = limits[1][1],extent=Map_plot)
	im2 = ax2.imshow(models_inv.B.transpose(), cmap=cmap[4], origin = origin, vmin = limits[4][0], vmax = limits[4][1],extent=Map_plot)
	im3 = ax3.imshow(models_inv.vlos.transpose(), cmap=cmap[5], origin = origin, vmin = limits[5][0], vmax = limits[5][1],extent=Map_plot)
	if "-plot_chi2" in sys.argv:
		im4 = ax4.imshow(chi2.tot.transpose(), cmap='Greys', origin = origin,extent=Map_plot)
	elif "-plot_fill" in sys.argv:
		im4 = ax4.imshow(models_inv.fill.transpose(),'Greys', origin = origin, vmin = 0, vmax = 1,extent=Map_plot)
	else:
		im4 = ax4.imshow(models_inv.gamma.transpose(), cmap=cmap[6], origin = origin, vmin = limits[6][0], vmax = limits[6][1],extent=Map_plot)
	if "-arc" in sys.argv and "-varrow" in sys.argv:
		l = np.min([Map_plot[1],Map_plot[3]])*0.025
		ax1.arrow(Map_plot[1]/2, Map_plot[3]/2, sign1*abs(dx)*n, sign2*abs(dy)*n, head_width=l, head_length=l, fc='black', ec='black')
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
		ax1.set_title(titles[1] + r" @ $\log\tau = $ " + str(logT))
		ax2.set_title(titles[4] + r" @ $\log\tau = $ " + str(logB))
		ax3.set_title(titles[5] + r" @ $\log\tau = $ " + str(logV))
		if "-plot_chi2" in sys.argv:
			ax4.set_title(titles[11])
		elif "-plot_fill" in sys.argv:
			ax4.set_title("Filling Factor")
		else:
			ax4.set_title(titles[6] + r" @ $\log\tau = $ " + str(logI))

	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.057 * frac, pad=0.04)
	cbar1.set_label(label = labels[1], loc = 'center', labelpad=15)
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.057 * frac, pad=0.04)
	cbar2.set_label(label = labels[4], loc = 'center', labelpad=15)
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.057 * frac, pad=0.04)
	cbar3.set_label(label = labels[5], loc = 'center', labelpad=15)
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.057 * frac, pad=0.04)
	if "-plot_chi2" in sys.argv:
		cbar4.set_label(label = labels[11], loc = 'center', labelpad=15)
	elif "-plot_fill" in sys.argv:
		cbar4.set_label(label = r"$\alpha$", loc = 'center', labelpad=15)	
	else:
		cbar4.set_label(label = labels[6], loc = 'center', labelpad=15)
	############
	# Set title position depending on the chosen plot and consider the flags hinode and gris
	if Type != "2":
		if title3 != "-1":
			if "-vertical" in sys.argv:
				xtitle1 = 0.41
			else:
				xtitle1 = 0.5
			if "-xtitle" in sys.argv:
				xtitle1 = float(sys.argv[sys.argv.index("-xtitle")+1])
			if title3 != '':
				fig.suptitle(title3, y=1.02, x=xtitle1)

	if Type == "2":
		if title4 != "-1":
			if "-vertical" in sys.argv:
				xtitle1 = 0.41
			else:
				xtitle1 = 0.5
			if "-xtitle" in sys.argv:
				xtitle1 = float(sys.argv[sys.argv.index("-xtitle")+1])
			if title4 != '':
				fig.suptitle(title4, y=1.02, x=xtitle1)
	

	plt.savefig(savepath + "inversion" + Type + add)



def _plot_stokes(stokes, stokes_inv, wave, Map, figsize, frac, units, title1,  title2, Map_plot, origin, sign1, sign2, dx, dy, n, savepath, add):
	'''
	Print the stokes vector (for internal use)
	'''

	waves = stokes.wave
	waves_inv = stokes_inv.wave

	# Determine indexes for the wavelength
	wave_ind1 = np.argmin(abs(waves-wave))
	wave_ind2 = np.argmin(abs(waves_inv-wave))
	wave = waves[wave_ind1]

	# Check whether the wavelength is in range
	wave = _check_range(stokes_inv.wave, wave)
	
	if "-waveV" in sys.argv:
		waveV = float(sys.argv[sys.argv.index("-waveV")+1])
		waveV = _check_range(stokes_inv.wave, waveV)
	else:
		waveV = wave
	waveV_ind1 = np.argmin(abs(stokes.wave-waveV))
	waveV_ind2 = np.argmin(abs(stokes_inv.wave-waveV))



	if "-waveQ" in sys.argv:
		waveQ = float(sys.argv[sys.argv.index("-waveQ")+1])
		waveQ = _check_range(stokes_inv.wave, waveQ)
	else:
		waveQ = wave
	waveQ_ind1 = np.argmin(abs(stokes.wave-waveQ))
	waveQ_ind2 = np.argmin(abs(stokes_inv.wave-waveQ))

	if "-waveU" in sys.argv:
		waveU = float(sys.argv[sys.argv.index("-waveU")+1])
		waveU = _check_range(stokes_inv.wave, waveU)
	else:
		waveU = wave
	waveU_ind1 = np.argmin(abs(stokes.wave-waveU))
	waveU_ind2 = np.argmin(abs(stokes_inv.wave-waveU))

	if not stokes._data_cut_map:
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
					[-np.max(np.abs(Q1[:,:,waveQ_ind1])), np.max(np.abs(Q1[:,:,waveQ_ind1]))],
					[-np.max(np.abs(U1[:,:,waveU_ind1])), np.max(np.abs(U1[:,:,waveU_ind1]))],
					[-np.max(np.abs(V1[:,:,waveV_ind1])), np.max(np.abs(V1[:,:,waveV_ind1]))]
		]
	limits_stokes2 = [  [np.min(I2[:,:,wave_ind2]), np.max(I2[:,:,wave_ind2])],
					[-np.max(np.abs(Q2[:,:,waveQ_ind2])), np.max(np.abs(Q2[:,:,waveQ_ind2]))],
					[-np.max(np.abs(U2[:,:,waveU_ind2])), np.max(np.abs(U2[:,:,waveU_ind2]))],
					[-np.max(np.abs(V2[:,:,waveV_ind2])), np.max(np.abs(V2[:,:,waveV_ind2]))]
		]


	if "-limitI" in sys.argv:
		temp = sys.argv[sys.argv.index("-limitI")+1].split(',')
		limits_stokes1[0] = [int(i) for i in temp ]
		limits_stokes2[0] = [int(i) for i in temp ]
	
		

	if "-vertical" in sys.argv:
		fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True,
			 gridspec_kw=dict(hspace=0), figsize=figsize)
		fig.subplots_adjust(hspace=0, wspace=0)
	else:
		if (I2.shape[1] - I2.shape[0]) >= -100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0]*2,figsize[1]*2/4],
													layout="compressed",
													)
			
		else:
			fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=figsize,
										layout="compressed",
									)

	############################
	# Plot the Stokes profiles #
	############################
	im1 = ax1.imshow(I1[:,:,wave_ind1]  .transpose(), origin=origin, vmin = limits_stokes1[0][0], vmax = limits_stokes1[0][1], extent=Map_plot)
	im2 = ax2.imshow(Q1[:,:,waveQ_ind1]  .transpose(), origin=origin, vmin = limits_stokes1[1][0], vmax = limits_stokes1[1][1], cmap = 'PuOr', extent=Map_plot)
	im3 = ax3.imshow(U1[:,:,waveU_ind1]  .transpose(), origin=origin, vmin = limits_stokes1[2][0], vmax = limits_stokes1[2][1], cmap = 'PuOr', extent=Map_plot)
	im4 = ax4.imshow(V1[:,:,waveV_ind1].transpose(), origin=origin, vmin = limits_stokes1[3][0], vmax = limits_stokes1[3][1], cmap = 'PuOr', extent=Map_plot)
	if "-arc" in sys.argv:
		l = np.min([Map_plot[1],Map_plot[3]])*0.025
		ax1.arrow(Map_plot[1]/2, Map_plot[3]/2, sign1*abs(dx)*n, sign2*abs(dy)*n, head_width=l, head_length=l, fc='black', ec='black')
	#####################
	#	Set labels	#
	#####################
	ax1.set_title(r'$\mathrm{I} / \mathrm{I}_c $ @' + "%.3f" % wave + r" \AA")
	ax2.set_title(r'$\mathrm{Q} / \mathrm{I}_c$ @' + "%.3f" % waveQ + r" \AA")
	ax3.set_title(r'$\mathrm{U} / \mathrm{I}_c$ @' + "%.3f" % waveU + r" \AA")
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
		ax1.set_xlabel(f"x [{units}]")
		ax2.set_xlabel(f"x [{units}]")
		ax3.set_xlabel(f"x [{units}]")
	ax4.set_xlabel(f"x [{units}]")
	ax1.set_ylabel(f"y [{units}]")
	ax2.set_ylabel(f"y [{units}]")
	ax3.set_ylabel(f"y [{units}]")
	ax4.set_ylabel(f"y [{units}]")

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
		if (I2.shape[1] - I2.shape[0]) >= -100:
			fig, ((ax1,ax2,ax3,ax4)) = plt.subplots(1,4,figsize=[figsize[0]*2,figsize[1]*2/4],
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
	im1 = ax1.imshow(I2[:,:,wave_ind2]  .transpose(), origin=origin, vmin = limits_stokes1[0][0], vmax = limits_stokes1[0][1], extent=Map_plot)
	im2 = ax2.imshow(Q2[:,:,waveQ_ind2]  .transpose(), origin=origin, vmin = limits_stokes1[1][0], vmax = limits_stokes1[1][1], cmap = 'PuOr', extent=Map_plot)
	im3 = ax3.imshow(U2[:,:,waveU_ind2]  .transpose(), origin=origin, vmin = limits_stokes1[2][0], vmax = limits_stokes1[2][1], cmap = 'PuOr', extent=Map_plot)
	im4 = ax4.imshow(V2[:,:,waveV_ind2].transpose(), origin=origin, vmin = limits_stokes1[3][0], vmax = limits_stokes1[3][1], cmap = 'PuOr', extent=Map_plot)

	
	if "-arc" in sys.argv:
		l = np.min([Map_plot[1],Map_plot[3]])*0.025
		ax1.arrow(Map_plot[1]/2, Map_plot[3]/2, sign1*abs(dx)*n, sign2*abs(dy)*n, head_width=l, head_length=l, fc='black', ec='black')

	#####################
	#	Set labels	#
	#####################
	ax1.set_title(r'$\mathrm{I} / \mathrm{I}_c $ @' + "%.3f" % wave + r" \AA")
	ax2.set_title(r'$\mathrm{Q} / \mathrm{I}_c$ @' + "%.3f" % waveQ + r" \AA")
	ax3.set_title(r'$\mathrm{U} / \mathrm{I}_c$ @' + "%.3f" % waveU + r" \AA")		
	ax4.set_title(r'$\mathrm{V} / \mathrm{I}_c$ @' + "%.3f" % waveV + r" \AA")		


	############
	# Colorbar #
	cbar1 = fig.colorbar(im1, ax=ax1, fraction=0.057 * frac, pad=0.04)
	cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.057 * frac, pad=0.04)
	cbar3 = fig.colorbar(im3, ax=ax3, fraction=0.057 * frac, pad=0.04)
	cbar4 = fig.colorbar(im4, ax=ax4, fraction=0.057 * frac, pad=0.04)
	############

	#####################
	#	Set labels	#
	#####################
	if "-vertical" not in sys.argv:
		ax1.set_xlabel(f"x [{units}]")
		ax2.set_xlabel(f"x [{units}]")
		ax3.set_xlabel(f"x [{units}]")
	ax4.set_xlabel(f"x [{units}]")
	ax1.set_ylabel(f"y [{units}]")
	ax2.set_ylabel(f"y [{units}]")
	ax3.set_ylabel(f"y [{units}]")
	ax4.set_ylabel(f"y [{units}]")

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

	plt.savefig(savepath + "stokes" + add)

def result(conf, wave, tau, Type = "", plot_stokes = True):
	"""
	Plots the result of the inversion

	Parameters
	----------
	conf : dict
		Dict. with all the information from the config
	wave : float
		Wavelength in A where the Stokes vector is plottet
	tau : float
		log tau value where the model is plotted
	Type : string, optional
		prefix for determining which model is used. Default: '' for Mode 1C, "1" for model 1 from mode 2C, "2" for model 2 from mode 2C
	plot_stokes : bool, optional
		Plot the stokes vector. Default: True

	Returns
	-------
	None

	Other Parameters
	----------------
	Additional parameters given as an argument when the script is executed.
	
	-data [str]
		Rel. path to the spectral veil corrected data if standard labelling is not used, optional.
	-stokes [str]
		Rel. path to the Stokes result if standard labelling is not used, optional.
	-model [str]
		Rel. path to the Models of the inversion if standard labelling is not used (Mode `1C`).
	-model1 [str]
		Rel. path to the Models 1 of the inversion if standard labelling is not used (Mode `2C`).
	-model2 [str]
		Rel. path to the Models 2 of the inversion if standard labelling is not used (Mode `2C`).
	-chi
		Rel. path to the chi2 file of the inversion if standard labelling is not used.
	-save [str], optional
		Additional save path. Default './'.
	-add [str]
		Additional text in filenames
	-label [str]
		Add label text
	-title1 [str]
		Title in Result Stokes plot
	-title2 [str]
		Title in Obs. Stokes plot
	-title3 [str]
		Title in Model (1) plot
	-title4 [str]
		Title in Model 2 plot with 4 plots (Mode `1C`)
	-xtitle [float]
		Changing the x position of the title
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
	-chi2
		Plot chi2
	-fill
		Plot the filling factor
	-plot_chi2
		Plot chi2 in the 4 subplots plot
	-plot_fill
		Plot filling factor in the 4 subplots plot
	-waveQ [float]
		Plot Stokes Q in another wavelength position in A
	-waveU [float]
		Plot Stokes U in another wavelength position in A
	-waveV [float]
		Plot Stokes V in another wavelength position in A
	-logT [float]
		Plot Temperature at this log tau.
	-logB [float]
		Plot Magnetic Field at this log tau.
	-logV [float]
		Plot LoS Velocity at this log tau.
	-logI [float]
		Plot Inclination at this log tau.
	-limitxy [float,float,float,float]
		Limit in the x and y plot as a list xmin,xmax,ymin,xmax
	-limitT [float,float]
		Set the limit for the colorbar in the temperature.
	-limitB [float,float]
		Set the limit for the colorbar in the magnetic field.
	-limitchi2 [float,float]
		Set the limit for the colorbar in chi2.
	-limitI [float,float]
		Set the limit for the colorbar in Stokes I.
	-kG
		Plot in kG (Kilogauss)
	-arc
		Print x and y axis in arcseconds
	-flipy
		Flip image along y axis as sometimes the orientation is wrong in GRIS with the location on the sun
	-f [float]
		Factor for the fontsizes
	-symv
		Symmetric limits for vlos
	-rot90
		Rotate the image 90 deg
	"""

	# Import library
	sir.mpl_library()
	
	# Check if path exists
	if not exists(conf['path']):
		print(f"[NOTE] Path {conf['path']} does not exist.")
		return

	#############################################################
	#			READ INPUT AND LOAD DATA					#
	#############################################################
	path = conf["path"]
	
	Map = conf['map']
	
	if "-data" not in sys.argv:
		stokes = p.read_profile(os.path.join(conf["path"],conf['cube']))
	else:
		filename = sys.argv[sys.argv.index("-data")+1]
		stokes = p.read_profile(filename)
	stokes = stokes.cut_to_map(conf["map"])
	stokes = stokes.cut_to_wave(conf["range_wave"])

	if "-stokes" not in sys.argv:
		stokes_inv = p.read_profile(os.path.join(path,conf['inv_out'] + d.end_stokes))
	else:
		filename = sys.argv[sys.argv.index("-stokes")+1]
		stokes_inv = p.read_profile(filename)
	

	if Type == "1":
		if "-model1" not in sys.argv:
			models_inv = m.read_model(os.path.join(path, conf['inv_out'] + d.end_models1))
		else:
			filename = sys.argv[sys.argv.index("-model1")+1]
			models_inv = m.read_model(filename)
	elif Type == "2":
		if "-model2" not in sys.argv:
			models_inv = m.read_model(os.path.join(path, conf['inv_out'] + d.end_models2))
		else:
			filename = sys.argv[sys.argv.index("-model2")+1]
			models_inv = m.read_model(filename)
	else:
		if "-model" not in sys.argv:
			models_inv = m.read_model(os.path.join(path, conf['inv_out'] + d.end_models))
		else:
			filename = sys.argv[sys.argv.index("-model")+1]
			models_inv = m.read_model(filename)

	if ("-chi2" in sys.argv or "-plot_chi2" in sys.argv):
		if "-chi" not in sys.argv:
			chi2 = c.read_chi2(os.path.join(path,conf['chi2']))
		else:
			filename = sys.argv[sys.argv.index("-chi")+1]
			chi2 = c.read_chi2(filename)
	else:
		chi2 = None

	# Rotate by 90 deg
	if "-rot90" in sys.argv:
		print("[INFO] Images are rotated by 90°!")
		models_inv.T = np.moveaxis(models_inv.T,(0,1,2),(1,0,2))
		models_inv.Pe = np.moveaxis(models_inv.Pe,(0,1,2),(1,0,2))
		models_inv.vmicro = np.moveaxis(models_inv.vmicro,(0,1,2),(1,0,2))
		models_inv.B = np.moveaxis(models_inv.B,(0,1,2),(1,0,2))
		models_inv.vlos = np.moveaxis(models_inv.vlos,(0,1,2),(1,0,2))
		models_inv.gamma = np.moveaxis(models_inv.gamma,(0,1,2),(1,0,2))
		models_inv.phi = np.moveaxis(models_inv.phi,(0,1,2),(1,0,2))
		models_inv.z = np.moveaxis(models_inv.z,(0,1,2),(1,0,2))
		models_inv.rho = np.moveaxis(models_inv.rho,(0,1,2),(1,0,2))
		models_inv.Pg = np.moveaxis(models_inv.Pg,(0,1,2),(1,0,2))
		models_inv.nx, models_inv.ny = models_inv.ny, models_inv.nx

		stokes.stki = np.moveaxis(stokes.stki,(0,1,2),(1,0,2))
		stokes.stkq = np.moveaxis(stokes.stkq,(0,1,2),(1,0,2))
		stokes.stku = np.moveaxis(stokes.stku,(0,1,2),(1,0,2))
		stokes.stkv = np.moveaxis(stokes.stkv,(0,1,2),(1,0,2))
		stokes.nx, stokes.ny = stokes.ny, stokes.nx

		stokes_inv.stki = np.moveaxis(stokes_inv.stki,(0,1,2),(1,0,2))
		stokes_inv.stkq = np.moveaxis(stokes_inv.stkq,(0,1,2),(1,0,2))
		stokes_inv.stku = np.moveaxis(stokes_inv.stku,(0,1,2),(1,0,2))
		stokes_inv.stkv = np.moveaxis(stokes_inv.stkv,(0,1,2),(1,0,2))
		stokes_inv.nx, stokes_inv.ny = stokes_inv.ny, stokes_inv.nx
		if ("-chi2" in sys.argv or "-plot_chi2" in sys.argv):
			chi2.stki = np.moveaxis(chi2.stki,(0,1,2),(1,0,2))
			chi2.stkq = np.moveaxis(chi2.stkq,(0,1,2),(1,0,2))
			chi2.stku = np.moveaxis(chi2.stku,(0,1,2),(1,0,2))
			chi2.stkv = np.moveaxis(chi2.stkv,(0,1,2),(1,0,2))
			chi2.tot = np.moveaxis(chi2.tot,(0,1,2),(1,0,2))
			chi2.nx, chi2.ny = chi2.ny, chi2.nx
		Map = [Map[2],Map[3],Map[0],Map[1]]


	# Cut data in x and y position	
	if "-limitxy" in sys.argv:
		limit_xy = np.array([int(i) for i in sys.argv[sys.argv.index("-limitxy")+1].split(",")], dtype=int)

		# Cut data to the new range:
		stokes_inv = stokes_inv.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])
		stokes = stokes.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])
		models_inv = models_inv.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])
		#errors_inv = errors_inv.cut_to_map([limit_xy[0]-Map[0],limit_xy[1]+1-(Map[1]+1), limit_xy[2]-Map[2],limit_xy[3]+1-(Map[3]+1)])

		# Save the new limits as the new Map
		Map = limit_xy


	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = path + "/" + sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath[:savepath.rfind('/')]):
			os.mkdir(savepath[:savepath.rfind('/')])

	# Additional text in output
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
	title4 = ''
	if '-title4' in sys.argv:
		title4 = sys.argv[sys.argv.index("-title4")+1]

	#########################
	#  PLOTTING STUFF		#
	#########################
	if (stokes_inv.ny - stokes_inv.nx) >= -100:
		alpha = stokes_inv.nx / 12
		figsize = [stokes_inv.nx / alpha , stokes_inv.ny/alpha-1]
		frac = figsize[1] / figsize[0]
		if "-f" in sys.argv:
			f = float(sys.argv[sys.argv.index("-f")+1])
		else:
			f = 1.2
	else: # for 2x2 plot
		alpha = stokes_inv.nx / 12
		figsize = [stokes_inv.nx / alpha, stokes_inv.ny/alpha]
		frac = figsize[1] / figsize[0] * 1.25
		if "-f" in sys.argv:
			f = float(sys.argv[sys.argv.index("-f")+1])
		else:
			f = 0.65
	
	import matplotlib as mpl
	mpl.rcParams["xtick.labelsize"] = 18*f
	mpl.rcParams["ytick.labelsize"] = 18*f
	mpl.rcParams["legend.fontsize"] = 16*f
	mpl.rcParams["legend.title_fontsize"] = 16*f
	mpl.rcParams["axes.titlesize"] = 18*f
	mpl.rcParams["axes.labelsize"] = 20*f
	mpl.rcParams["figure.titlesize"] = 24*f

	####################################
	#	Plot settings for arcsecond	#
	####################################
	dx = dy = sign1 = sign2 = n = None
	if "-arc" in sys.argv:
		infos = dict(np.genfromtxt(d.header_infos,dtype='str', delimiter="="), dtype=str)
		if conf['instrument'] == 'GRIS':
			x = float(infos['CRVAL1'])
			y = float(infos['CRVAL2'])
			dx = float(infos['CDELT1'])
			dy = float(infos['CDELT2'])	
		elif conf['instrument'] == 'Hinode':
			dy = float(infos['CDELT2'])  # Delta y of slit
			dx = float(infos['XSCALE'])
			x = float(infos['XCEN'])
			y = float(infos['YCEN'])
			
		Map_plot = [0,stokes_inv.nx*abs(dx),0,stokes_inv.ny*abs(dy)]
		if x > 0 and y > 0:
			sign1 = -1
			sign2 = -1
		elif x > 0 and y < 0:
			sign1 = -1
			sign2 = +1
		elif x < 0 and y > 0:
			sign1 = +1
			sign2 = -1
		elif x < 0 and y < 0:
			sign1 = +1
			sign2 = +1
		
		# Account for the angle between (x,y) and the center
		alpha = np.arctan(abs((y+Map_plot[3]/2)/(x+Map_plot[1]/2))) # Angle at the center of the image
		sign1 = sign1*np.cos(alpha) # abs because direction is handled before
		sign2 = sign2*np.sin(alpha) # abs because direction is handled before

		n = np.min([stokes.nx,stokes.ny])*0.2
	else:
		Map_plot = Map
	
	origin = d.origin

	if "-flipy" in sys.argv:
		print("Image flipped in y direction")
		if origin == "lower":
			origin = "upper"
		elif origin == "upper":
			origin = "lower"
		
	if "-arc" in sys.argv:
		units = 'Arcsec'
	else:
		units = 'Pixels'

	if plot_stokes:
		_plot_stokes(stokes, stokes_inv, wave, Map, figsize, frac, units, title1,  title2, Map_plot, origin, sign1, sign2, dx, dy, n, savepath, add)

	###################################################
	#			Plot physical parameters			#
	###################################################
	_plot_model(models_inv, tau, figsize, frac, units, title3, title4, savepath, add, chi2, Map_plot, origin, sign1, sign2, n, dx, dy, Type)

	if "-chi2" in sys.argv:
		plot_chi2(figsize, frac, chi2, Map_plot, units, savepath, add, origin)
# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		_help()
	conf = sir.read_config(sys.argv[1])

	if conf['mode'] == '1C':
		result(conf, float(sys.argv[2]), float(sys.argv[3]),"",True)

	elif conf['mode'] == '2C':
		result(conf, float(sys.argv[2]), float(sys.argv[3]), "1")
		result(conf, float(sys.argv[2]), float(sys.argv[3]), Type = "2" , plot_stokes = False)

	else:
		print(f"[result] Mode '{conf['mode']}' not known or undefined!")





