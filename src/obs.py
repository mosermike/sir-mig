"""
Functions related to observations and to conversions between SIR format and numpy files
"""
import numpy as np

def read_profile(profile, grid, line_file):
	"""
	Reads a profile from a SIR .per file

	Parameters
	----------
	profile : str
		String with the profile name
	grid : str
		Grid file name
	line_file : str
		Filename of the line file

	Returns
	-------
	out : numpy.array
		Array with the wavelengths
	out : numpy.array
		Array with Stokes I
	out : numpy.array
		Array with Stokes Q
	out : numpy.array
		Array with Stokes U
	out : numpy.array
		Array with Stokes V

	"""
	import sir
	# Read grid and line file to transform the relative llambdas
	grid = sir.read_grid(grid)
	Lines = grid['Line']

	#Line file
	line = sir.read_line(line_file)
	wave = line['wavelength']
	line = line['Line']

	
	# Read the profile
	lls = np.array([])
	Is = np.array([])
	Qs = np.array([])
	Us = np.array([])
	Vs = np.array([])
	for i in Lines:
		ll, I, Q, U, V = sir.read_profile(profile, num = int(i[0]))
		lls = np.append(lls, np.array(ll) + wave[np.where(line == int(i[0]))]) # in abs. wavelength
		Is = np.append(Is, I)
		Qs = np.append(Qs, Q)
		Us = np.append(Us,U)
		Vs = np.append(Vs, V)

	return lls, Is, Qs, Us, Vs


##############################################################################################








