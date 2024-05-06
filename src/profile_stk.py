"""
Class Profile with all the tools to read and write the Stokes Profiles

"""

import numpy as np 
import os
import sys
from scipy.io import FortranFile


class Profile:
	"""
	Class containing the models of a simulation

	Variables are:
		- wave
		- indx (= line numbers in mode 'MC')
		- stki
		- stkq
		- stku
		- stkv


	There are several functions:
		- read: Read a numpy file containing models and stores it into the class variables
		- read_results: Reads the results from the inversion of SIR
		- write_profile: Writes a Model to a SIR readable file
		- write: Saves the models as a numpy file
		- set_limit: Cuts the data to a specific log_tau value (used for plotting)


	"""
	def __init__(self, nx=None, ny = None, nw=0):
		"""
		Initialisation of the class with the Profiles
		
		Parameter
		---------
		nx : int
			Integer of pixels in x direction
		ny :  int
			Integer of pixels in y direction
		nw : int
			Number of wavelength points
		"""
		filename = None
		# Initialize with two integers
		self.nx = nx	# Points in x
		self.ny = ny	# Points in y
		self.nw = nw
		self.ns = 4
		
		self.load = False		# Determines whether data is already loaded or not
		
		self.indx = np.zeros(shape=(self.nw))
		self.wave = np.zeros(shape=(self.nw))
		self.stki = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float32)
		self.stkq = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float32)
		self.stku = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float32)
		self.stkv = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float32)

		self.data_cut = False

	def __read_grid(self, filename):
		"""
		Reads the grid file
		
		Parameter
		---------
		filename : string
			File to be read
		
		Return
		-------
		dict : Dictionary
			Dict. with 'Line', 'min', 'step' and 'max' in it
		"""
		# Open the file and read lines
		with open(filename) as f:
			strings = f.readlines()

		# Remove last line if it contains no information
		if ":" not in strings[-1]:
			strings = strings[:-1]

		# Remove leading spaces
		for i in range(len(strings)):
			while(strings[i][0] == ' '):
				strings[i] = strings[i][1:]
		
		# Create arrays
		Line	  = []
		ll_min	= np.empty(len(strings))
		ll_step    = np.empty(len(strings))
		ll_max	= np.empty(len(strings))

		
		# Remove multiple spaces and split 
		strings = [i.replace(':', ': ')   for i in strings] # Add space to make sure it is split correctly
		while any('  ' in x for x in strings):
			strings = [i.replace('  ', ' ')   for i in strings]
		strings = [i.split(' ') for i in strings]

		# Fill arrays with data	
		for i in range(len(strings)):
			Line.append(strings[i][0][:-1].split(','))
			ll_min[i]  = float(strings[i][1].replace(",",""))
			ll_step[i] = float(strings[i][2].replace(",",""))
			ll_max[i]  = float(strings[i][3].replace(",",""))

		Dict = {
				'Line'		: Line,
				'min'		: ll_min,
				'step'		: ll_step,
				'max'		: ll_max,
			}
		
		return Dict


	def __read_profile_sir(self, filename):
		"""
		Reads the first LINE data from a profile computed by SIR
		
		Parameter
		---------
		filename : string
			String containing the path of the file
		
		Return
		------
		ll : numpy.array
			Wavelengths in A
		I : numpy.array
			Stokes I
		Q : numpy.array
			Stokes Q
		U : numpy.array
			Stokes U
		V : numpy.array 
			Stokes V
		"""
		data = np.loadtxt(filename).transpose()
		line = data[0].astype(int)
		ll = data[1].astype(np.float32)
		I  = data[2].astype(np.float32)
		Q  = data[3].astype(np.float32)
		U  = data[4].astype(np.float32)
		V  = data[5].astype(np.float32)

		return np.array(ll), np.array(I), np.array(Q), np.array(U), np.array(V)

	def __read_profile_sir_mc(self, filename):
		"""
		Reads the first LINE data from a profile computed by SIR
		
		Parameter
		---------
		filename : string
			String containing the path of the file
		
		Return
		------
		ll : numpy.array
			Wavelengths in A
		I : numpy.array
			Stokes I
		Q : numpy.array
			Stokes Q
		U : numpy.array
			Stokes U
		V : numpy.array 
			Stokes V
		"""
		if not os.path.exists(filename):
			print(f"[ERROR] File {filename} does not exist.")

		data = np.loadtxt(filename).transpose()
		line = data[0].astype(int)
		ll = data[1].astype(np.float32)
		I  = data[2].astype(np.float32)
		Q  = data[3].astype(np.float32)
		U  = data[4].astype(np.float32)
		V  = data[5].astype(np.float32)

		return line, np.array(ll), np.array(I), np.array(Q), np.array(U), np.array(V)


	def cut_to_map(self, Map):
		"""
		Cut the data to a map [xmin, xmax, ymin, ymax]

		Parameters
		----------
		Map : list
			List with the ranges in pixel in x and y direction
		"""
		
		self.nx = Map[1]-Map[0]+1
		self.ny = Map[3]-Map[2]+1

		self.stki = self.stki[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.stkq = self.stkq[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.stku = self.stku[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.stkv = self.stkv[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		
		
		return self

	def cut_to_wave(self, range_wave):
		"""
		Cut the data to the range in wavelengths

		Parameters
		----------
		range_wave : list
			List with the ranges from the config file
		"""
		# Number of wavelengths
		temp = range_wave[:,2].astype(int)
		nws = [0, temp[0]]
		for i in range(1,len(temp)):
			nws.append(temp[i]+temp[i-1])
		
		# Initialize new arrays
		lwave = np.zeros(shape=(nws[-1]), dtype=np.float32)
		lstki = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float32)
		lstkq = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float32)
		lstku = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float32)
		lstkv = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float32)
		
		ind = [np.argmin(np.abs(self.wave-range_wave[i][0])) for i in range(len(range_wave))]
		for i in range(len(range_wave)):
			lwave[nws[i]:nws[i+1]] = self.wave[ind[i]:ind[i]+int(range_wave[i][2])]
			lstki[:,:,nws[i]:nws[i+1]] = self.stki[:,:,ind[i]:ind[i]+int(range_wave[i][2])]
			lstkq[:,:,nws[i]:nws[i+1]] = self.stkq[:,:,ind[i]:ind[i]+int(range_wave[i][2])]
			lstku[:,:,nws[i]:nws[i+1]] = self.stku[:,:,ind[i]:ind[i]+int(range_wave[i][2])]
			lstkv[:,:,nws[i]:nws[i+1]] = self.stkv[:,:,ind[i]:ind[i]+int(range_wave[i][2])]
		
		self.wave = lwave
		self.stki = lstki
		self.stkq = lstkq
		self.stku = lstku
		self.stkv = lstkv

		self.nw = self.wave.shape[0]
		self.data_cut = True
		
		return self
	

	def read(self, fname, fmt_type=np.float32):
		
		f = FortranFile(fname, 'r')
		first_rec = f.read_record(dtype=fmt_type)

		posv = first_rec[0]
		negv = first_rec[1]
		fid = first_rec[2]
		nrec = first_rec[3]
		ndims = first_rec[4]

		medv = (posv + negv) / 2.
		verp = posv - medv
		vern = negv - medv
		vers = (posv - negv) / 2.

		if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-160904.) > 1.e-3) ):
			print('Something is wrong with the file %s' % fname)
			print('\tIs it a 3D profile file?')
			return np.nan

		if (np.abs(vers-3.) < 1.e-3):
			#VERSION 3
			dims = first_rec[5:]
			nx, ny, nw, ns = np.int16(dims)
			self.nx = nx
			self.ny = ny
			self.nw = nw
			self.ns = ns

			self.indx = f.read_record(dtype=fmt_type).astype(np.float64)

			self.wave = f.read_record(dtype=fmt_type).astype(np.float64)

			data = f.read_record(dtype=fmt_type)
			data = data.reshape(nx,ny,nw,ns) * 1.

			self.stki = data[:, :, :, 0].astype(np.float64)
			self.stkq = data[:, :, :, 1].astype(np.float64)
			self.stku = data[:, :, :, 2].astype(np.float64)
			self.stkv = data[:, :, :, 3].astype(np.float64)
		
				
		self.load = True

		return self

	def read_results(self, task, filename, path, nx, ny):
		"""
		Reads all the errors from the inversion
		
		Parameter
		---------
		task : dict
			Dictionary with all the task folders, x and y positions
		filename : str
			Filename of the file in each task folder
		path : str
			Path where all the files are
		nx : int
			how many results are read in x
		ny : int
			how many results are read in y
		

		Return
		------
		class with the parameters
		"""	
		# Set the dimensions
		self.nx = nx
		self.ny = ny
		
		# Read the profile
		file = f"{os.path.join(path,task['folders'][0])}/{filename}"
		ll, _, _, _, _ = self.__read_profile_sir(file)
		self.nw = len(ll)
		self.stki = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		self.stkq = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		self.stku = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		self.stkv = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		
		# Read all the files
		for i in range(len(task['x'])):
			file = f"{os.path.join(path,task['folders'][i])}/{filename}"
			ll, I, Q, U, V = self.__read_profile_sir(file)
			if( (i == 0) and (self.wave == 0).all()):
				self.wave = ll
			x, y = task['x'][i]-task['x'][0], task['y'][i]-task['y'][0]
			self.stki[x, y, :] = I
			self.stkq[x, y, :] = Q
			self.stku[x, y, :] = U
			self.stkv[x, y, :] = V
			
		self.indx = np.zeros(self.wave.shape)
		self.load = True

		return self


	def read_results_MC(self, path, tasks, filename):
		"""
		Reads all the profiles of the synthesis or inversion

		config : dict
			Config information
		tasks : dict
			Dictionary with the folder names
		Type : string, optional
			Indicating if synthesis or inversion. Default: 'syn'

		"""
			
		# Make a check
		if not os.path.exists(os.path.join(path,tasks['folders'][0]) + '/' + filename):
			print('[read_profiles] The profiles do not exist. Make sure, that sir is executed correctly and fortran is installed.')
			sys.exit()

		line, _, _, _, _, _ = self.__read_profile_sir_mc(os.path.join(path,tasks['folders'][0]) + '/' + filename)
		self.nw = len(line)
		self.nx = len(tasks['folders'])
		self.ny = 1
		self.stki = np.zeros(shape=(self.nx,1,self.nw))
		self.stkq = np.zeros(shape=(self.nx,1,self.nw))
		self.stku = np.zeros(shape=(self.nx,1,self.nw))
		self.stkv = np.zeros(shape=(self.nx,1,self.nw))

		for n in range(self.nx):
			line, ll, I, Q, U, V = self.__read_profile_sir_mc(os.path.join(path,tasks['folders'][n]) + '/' + filename)
			if n == 0:
				self.indx = line
				self.wave = ll
			self.stki[n,0,:] = I
			self.stkq[n,0,:] = Q
			self.stku[n,0,:] = U
			self.stkv[n,0,:] = V

		del line
		del ll
		del I
		del Q
		del U
		del V
		
		self.load = True

		return self


	def set_dim(self, nx, ny, nw):
		"""
		Sets the dimensions if no data is loaded yet

		Parameter
		---------
		nx : int
			Number of Models in x
		ny : int
			Number of Models in y
		nw : int
			Number of wavelength points
		
		"""
		if self.load:
			print("[set_dim] Data is already loaded and therefore dimensions cannot be set.")
			return self

		self.indx = np.zeros(shape=(nw), dtype=np.float64)
		self.wave = np.zeros(shape=(nw), dtype=np.float64)
		self.stki = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stkq = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stku = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stkv = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		
		self.nx = nx
		self.ny = ny
		self.nw = nw
		self.ns = 4

		return self

	def veil_correction(self, nu, Ic = 1.):
		"""
		Correct the spectrum for the spectral veil by simly inverting
		the convoluted equation.

		Parameters
		----------
		nu : float
			Optimized fraction of spectral veil
		Ic : float, optional
			Continuum intensity of the FTS. Default is "1.0".
		"""
		self.stki =  (self.stki - nu*Ic) / (1 - nu)
		
		return self

	def write(self, fname, fmt_type=np.float32):
		"""
		Write data into a binary fortran file

		Parameter
		---------
		fname : str
			File name 
		fmt_type : type
			type which is used to save it => numpy.float64 used
		"""
		if (fmt_type != np.float32):
			print('Not implemented yet!')
			return np.nan

		if not self.load:
			print(f"[write] It seems that data was never read or saved. {fname} may contain no data")


		f = FortranFile(fname, 'w')
		
		# FIRST, WE WRITE A FIRST RECORD WITH THE MANDATORY DATA:
		towrite = np.zeros(9, dtype=fmt_type)
		
		v = 3

		towrite[0] = 3000 + v
		towrite[1] = 3000 - v
		towrite[2] = 160904. * 1.
		towrite[3] = 4. #NREC
		towrite[4] = 4. # Number of dimensions
		towrite[5] = self.nx * 1.
		towrite[6] = self.ny * 1.
		towrite[7] = self.nw * 1.
		towrite[8] = 4. # Number of Stokes parameter
		
		# Write header
		f.write_record(np.float32(towrite))
		
		# Write indx
		f.write_record(np.float32(self.indx))

		# Write the wavelenghts
		f.write_record(np.float32(self.wave))

		array = np.zeros(shape=(self.nx, self.ny, self.nw, 4))


		array[:,:,:,0] = self.stki
		array[:,:,:,1] = self.stkq
		array[:,:,:,2] = self.stku
		array[:,:,:,3] = self.stkv
		

		f.write_record(np.float32(array.flatten()))

		del array

		f.close()

		return self


	def write_profile(self, filename, x, y, Grid):
		"""
		Writes data to profiles as described in the config file
		Note to call cut_to_wave, otherwise it writes the wrong profiles!
		
		Parameter
		---------
		filename : string
			String containing the output path of the profile
		x : int
			Position in x
		y : int
			Position in y
		Grid : string
			Grid file
		"""
		if not self.data_cut:
			print("[WARN] The data was not cut. Consider running Profile.cut_to_wave! Make sure that only data used in the inversion is loaded!")

		# Load data from config
		grid = self.__read_grid(Grid) # Grid file
				
		# Read the grid
		Line      = grid['Line']
		Line_min  = grid['min']
		Line_step = grid['step']
		Line_max  = grid['max']
		
		add = np.copy(Line_step) # in case the Line_step does not need to be added
		
		# Create the first column => number of the line
		num = np.empty(0)
		checks = []
		for i in range(len(Line)):
			num = np.append(num,np.ones(int(np.ceil((Line_max[i]-Line_min[i]+add[i])/Line_step[i])))*int(Line[i][0]))
			checks.append(int(np.ceil((Line_max[i]-Line_min[i]+add[i])/Line_step[i])))

		# Check if the data has the same format as the cut data (in case because if precision it is one value more)
		if np.sum(checks) != self.wave.shape[0]:
			add = np.zeros(add.shape)
			num = np.empty(0)
			checks = []
			for i in range(len(Line)):
				num = np.append(num,np.ones(int(np.ceil((Line_max[i]-Line_min[i]+add[i])/Line_step[i])))*int(Line[i][0]))
				checks.append(int(np.ceil((Line_max[i]-Line_min[i]+add[i])/Line_step[i])))
		
		
		# Create the second column => wavelength grid
		ll = np.empty(0)
		
		for i in range(len(Line)):
			if np.arange(Line_min[i],Line_max[i]+add[i],Line_step[i]).shape[0] != checks[i]:
				print(f"[WARN] The profiles might we wrong written {np.arange(Line_min[i],Line_max[i]+add[i],Line_step[i]).shape[0]} vs {checks[i]}")
			ll = np.append(ll,np.arange(Line_min[i],Line_max[i]+add[i],Line_step[i]))

		# Create arrays of the Stokes vector in the given wavelength range
		# Cut data as needed
		I = self.stki[x,y]
		Q = self.stkq[x,y]
		U = self.stku[x,y]
		V = self.stkv[x,y]

		# Save data
		f = open(filename, 'w')
		for i in range(len(num)):
			f.write(f" {int(num[i]):>2} {ll[i]:>10.4f} {I[i]:>14.7E} {Q[i]:>14.7E} {U[i]:>14.7E} {V[i]:>14.7E}\n")
		f.close()


	def write_profile_mc(self, filename, x):
		"""
		Writes data to profiles as described in the config file
		
		Parameter
		---------
		filename : string
			String containing the output path of the profile
		x : int
			Position in x
		y : int
			Position in y
		Grid : string
			Grid file
		"""

		# Save data
		f = open(filename, 'w')
		for i in range(self.nw):
			f.write(f" {int(self.indx[i]):>2} {self.wave[i]:>10.4f} {self.stki[x,0,i]:>14.7E} {self.stkq[x,0,i]:>14.7E} {self.stku[x,0,i]:>14.7E} {self.stkv[x,0,i]:>14.7E}\n")
		f.close()


def read_profile(file):
	"""
	Reads a profile and returns a class

	Parameters
	----------
	file : str
		String with the filename

	Returns
	-------
	class Profile

	"""
	pro = Profile(0,0,0)
	pro.read(file)

	return pro

