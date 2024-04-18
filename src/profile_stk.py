"""
Class Profile

"""

import numpy as np 
import os
import sys

def read_grid(filename):
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


def read_profile(filename):
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
	ll = data[1].astype(np.float64)
	I  = data[2].astype(np.float64)
	Q  = data[3].astype(np.float64)
	U  = data[4].astype(np.float64)
	V  = data[5].astype(np.float64)

	return np.array(ll), np.array(I), np.array(Q), np.array(U), np.array(V)

def read_profile_mc(filename):
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
	ll = data[1].astype(np.float64)
	I  = data[2].astype(np.float64)
	Q  = data[3].astype(np.float64)
	U  = data[4].astype(np.float64)
	V  = data[5].astype(np.float64)

	return line, np.array(ll/1000), np.array(I), np.array(Q), np.array(U), np.array(V)

class Profile:
	"""
	Class containing the models of a simulation

	Variables are:
		- wave
		- stki
		- stkq
		- stku
		- stkv


	There are several functions:
		- read: Read a numpy file containing models and stores it into the class variables
		- read_results: Reads the results from the inversion of SIR
		- write: Writes a Model to a SIR readable file
		- save: Saves the models as a numpy file
		- set_limit: Cuts the data to a specific log_tau value (used for plotting)


	"""
	def __init__(self, arg1=None, arg2 = None, nw=0):
		"""
		Initialisation of the class with the Profiles
		
		Parameter
		---------
		filename : str or int
			File name of the profile to be loaded. Default: None (= no reading)
			Integer of pixels in x direction
		filename_wave : str or int
			File name of the wavelength to be loaded. Default: None (= no reading)
			Integer of pixels in y direction
		nw : int
			Number of wavelength points (only used if arg1 and arg2 are integers)
		"""
		filename = None
		filename_wave = None
		# Initialize with two integers
		if isinstance(arg1, int) and isinstance(arg2, int):
			self.nx = arg1	# Points in x
			self.ny = arg2	# Points in y
		# Initialize with two strings
		elif isinstance(arg1, str) and isinstance(arg2, str):
			filename = arg1
			filename_wave = arg2
		elif isinstance(arg1, str) and arg2 is None:
			filename = arg1
			filename_wave = None
		else:
			self.nx = 0
			self.ny = 0
		
		self.full = False		# Determines whether also z, pg, and rho exists
		self.load = False		# Determines whether data is already loaded or not
		self.nw = nw			# Points in lambda

		if filename_wave is not None:
			self.wave = np.load(filename_wave)
		else:
			self.wave = np.zeros(shape=(self.nw))

		if filename is not None:
			if filename[-4:] == ".per":
				ll, I, Q, U, V = read_profile(filename)
				self.stki = np.zeros(shape=(1,1,len(ll)))
				self.stkq = np.zeros(shape=(1,1,len(ll)))
				self.stku = np.zeros(shape=(1,1,len(ll)))
				self.stkv = np.zeros(shape=(1,1,len(ll)))
				self.stki[0,0] = I
				self.stkq[0,0] = Q
				self.stku[0,0] = U
				self.stkv[0,0] = V
				self.nx = 1
				self.ny = 1
				self.nw = len(I)
			else:
				stk = np.load(filename).astype(np.float64)
				self.stki = stk[:, :, 0, :]
				self.stkq = stk[:, :, 1, :]
				self.stku = stk[:, :, 2, :]
				self.stkv = stk[:, :, 3, :]
				self.nx = self.stki.shape[0]
				self.ny = self.stki.shape[1]
				self.nw = self.stki.shape[2]
			
			self.load = True

		else:
			self.stki = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float64)
			self.stkq = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float64)
			self.stku = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float64)
			self.stkv = np.zeros(shape=(self.nx,self.ny,self.nw), dtype=np.float64)

		self.data_cut = False
	
	



	def cut_to_wave(self, range_wave_pix):
		"""
		Cut the data to the range in wavelengths

		Parameters
		----------
		range_wave_pix : list
			List with the ranges in pixel
		"""
		nws = [0]
		
		# Determine number of points and then where the data is added
		for i in range(len(range_wave_pix)):
			nws.append(range_wave_pix[i][1]-range_wave_pix[i][0]+1+nws[i])
		
		# Initialize new arrays
		lwave = np.zeros(shape=(nws[-1]), dtype=np.float64)
		lstki = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float64)
		lstkq = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float64)
		lstku = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float64)
		lstkv = np.zeros(shape=(self.nx, self.ny,nws[-1]), dtype=np.float64)
		
		for i in range(len(range_wave_pix)):
			lwave[nws[i]:nws[i+1]] = self.wave[range_wave_pix[i][0]:range_wave_pix[i][1]+1]
			lstki[:,:,nws[i]:nws[i+1]] = self.stki[:,:,range_wave_pix[i][0]:range_wave_pix[i][1]+1]
			lstkq[:,:,nws[i]:nws[i+1]] = self.stkq[:,:,range_wave_pix[i][0]:range_wave_pix[i][1]+1]
			lstku[:,:,nws[i]:nws[i+1]] = self.stku[:,:,range_wave_pix[i][0]:range_wave_pix[i][1]+1]
			lstkv[:,:,nws[i]:nws[i+1]] = self.stkv[:,:,range_wave_pix[i][0]:range_wave_pix[i][1]+1]

		self.wave = lwave
		self.stki = lstki
		self.stkq = lstkq
		self.stku = lstku
		self.stkv = lstkv

		
		self.data_cut = True
		
		return self

	def read(self, filename, filename_wave):
		"""
		Reads a numpy file and stores the parameters
		
		Parameter
		---------
		filename : str
			File name of the model to be loaded
		filename_wave : str
			Filename with the wavelengths
		
		"""
		stk = np.load(filename).astype(np.float64)

		self.stki = stk[:, :, 0, :]
		self.stkq = stk[:, :, 1, :]
		self.stku = stk[:, :, 2, :]
		self.stkv = stk[:, :, 3, :]
		
		self.wave = np.load(filename_wave)
		
		self.load = True
		self.nx = self.stki.shape[0]
		self.ny = self.stki.shape[1]
		self.nw = self.stki.shape[2]

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
		ll, _, _, _, _ = read_profile(file)
		self.nw = len(ll)
		self.stki = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		self.stkq = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		self.stku = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		self.stkv = np.zeros(shape=(nx, ny, len(ll)), dtype=np.float64)
		
		# Read all the files
		for i in range(len(task['x'])):
			file = f"{os.path.join(path,task['folders'][i])}/{filename}"
			ll, I, Q, U, V = read_profile(file)
			if( (i == 0) and (self.wave == 0).all()):
				self.wave = ll
			x, y = task['x'][i]-task['x'][0], task['y'][i]-task['y'][0]
			self.stki[x, y, :] = I
			self.stkq[x, y, :] = Q
			self.stku[x, y, :] = U
			self.stkv[x, y, :] = V
			

		self.load = True

		return self

	def save(self, filename, filename_wave):
		"""
		Save the data in the class as a numpy file

		Parameters
		----------
		filename : str
			Filename of the numpy file
		filename_wave : str
			Filename with the wavelengths

		"""
		if not self.load:
			print(f"[save] It seems that data was never read or saved. {filename} may contain no data")
			
		stk = np.stack([self.stki, self.stkq, self.stku, self.stkv], axis=2, dtype=np.float64)
		np.save(filename_wave,self.wave)
		np.save(filename, stk)

		del stk  # make the memory free

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

		self.wave = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stki = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stkq = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stku = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		self.stkv = np.zeros(shape=(nx, ny, nw), dtype=np.float64)
		
		self.nx = nx
		self.ny = ny
		self.nw = nw

		return self


	def write(self, filename, x, y, Grid):
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
		grid = read_grid(Grid) # Grid file
				
		# Read the grid
		Line      = grid['Line']
		Line_min  = grid['min']
		Line_step = grid['step']
		Line_max  = grid['max']

		# Create the first column => number of the line
		num = np.empty(0)
		checks = []
		for i in range(len(Line)):
			num = np.append(num,np.ones(int(np.ceil((Line_max[i]-Line_min[i]+Line_step[i])/Line_step[i])))*int(Line[i][0]))
			checks.append(int(np.ceil((Line_max[i]-Line_min[i]+Line_step[i])/Line_step[i])))

		# Create the second column => wavelength grid
		ll = np.empty(0)
		
		for i in range(len(Line)):
			if np.arange(Line_min[i],Line_max[i]+Line_step[i],Line_step[i]).shape[0] != checks[i]:
				print(f"[WARN] The profiles might we wrong written {np.arange(Line_min[i],Line_max[i]+Line_step[i],Line_step[i]).shape[0]} vs {checks[i]}")
			ll = np.append(ll,np.arange(Line_min[i],Line_max[i]+Line_step[i],Line_step[i]))

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



class Profile_MC:
	"""
	Class containing the models of a simulation

	Variables are:
		- line
		- wave
		- stki
		- stkq
		- stku
		- stkv


	There are several functions:
		- read: Read a numpy file containing models and stores it into the class variables
		- read_results: Reads the results from the inversion of SIR
		- write: Writes a Model to a SIR readable file
		- save: Saves the models as a numpy file
		- set_limit: Cuts the data to a specific log_tau value (used for plotting)


	"""
	def __init__(self, arg1=None, arg2 = 0):
		"""
		Initialisation of the class with the Profiles
		
		Parameter
		---------
		filename : str or int
			File name of the profile to be loaded. Default: None (= no reading)
			Integer of pixels in x direction
		filename_wave : str or int
			File name of the wavelength to be loaded. Default: None (= no reading)
			Integer of pixels in y direction
		nw : int
			Number of wavelength points (only used if arg1 and arg2 are integers)
		"""
		filename = None
		# Initialize with two integers
		if isinstance(arg1, int):
			self.num = arg1	# Points in x
			self.nw = arg2	# Points in y
		# Initialize with two strings
		elif isinstance(arg1, str):
			filename = arg1
			self.nw = 0
		else:
			self.nx = 0
			self.nw = 0
		
		self.load = False		# Determines whether data is already loaded or not
		self.npars = 6

		if filename is not None:
			stk = np.load(filename).astype(np.float64)
			self.line = stk[:, 0, :]
			self.wave = stk[:, 1, :]
			self.stki = stk[:, 2, :]
			self.stkq = stk[:, 3, :]
			self.stku = stk[:, 4, :]
			self.stkv = stk[:, 5, :]
		
			self.load = True
			self.num = stk.shape[0]
			self.npars = stk.shape[1]
			self.nw = stk.shape[2]
		else:
			self.line = np.zeros(shape=(self.num,self.nw), dtype=np.float64)
			self.wave = np.zeros(shape=(self.num,self.nw), dtype=np.float64)
			self.stki = np.zeros(shape=(self.num,self.nw), dtype=np.float64)
			self.stkq = np.zeros(shape=(self.num,self.nw), dtype=np.float64)
			self.stku = np.zeros(shape=(self.num,self.nw), dtype=np.float64)
			self.stkv = np.zeros(shape=(self.num,self.nw), dtype=np.float64)
			

	def read(self, filename):
		"""
		Reads a numpy file and stores the parameters
		
		Parameter
		---------
		filename : str
			File name of the model to be loaded
		filename_wave : str
			Filename with the wavelengths
		
		"""
		stk = np.load(filename).astype(np.float64)

		self.line = stk[:,0,:]
		self.wave = stk[:,1,:]
		self.stki = stk[:, 2, :]
		self.stkq = stk[:, 3, :]
		self.stku = stk[:, 4, :]
		self.stkv = stk[:, 5, :]
		
		self.load = True
		self.num = self.stki.shape[0]
		self.npars = stk.shape[1]
		self.nw = self.stki.shape[1]

		return self
	

	def read_results(self, path, tasks, filename):
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

		line, _, _, _, _, _ = read_profile_mc(os.path.join(path,tasks['folders'][0]) + '/' + filename)
		self.nw = len(line)
		self.num = len(tasks['folders'])
		stokes = np.zeros(shape=(self.num,self.nw))
		self.line = np.zeros(shape=(self.num,self.nw))
		self.wave = np.zeros(shape=(self.num,self.nw))
		self.stki = np.zeros(shape=(self.num,self.nw))
		self.stkq = np.zeros(shape=(self.num,self.nw))
		self.stku = np.zeros(shape=(self.num,self.nw))
		self.stkv = np.zeros(shape=(self.num,self.nw))

		for n in range(self.num):
			line, ll, I, Q, U, V = read_profile_mc(os.path.join(path,tasks['folders'][n]) + '/' + filename)
			self.line[n,:] = line
			self.wave[n,:] = ll
			self.stki[n,:] = I
			self.stkq[n,:] = Q
			self.stku[n,:] = U
			self.stkv[n,:] = V

		del line
		del ll
		del I
		del Q
		del U
		del V
		
		self.load = True

		return stokes

	def save(self, filename):
		"""
		Save the data in the class as a numpy file

		Parameters
		----------
		filename : str
			Filename of the numpy file
		filename_wave : str
			Filename with the wavelengths

		"""
		if not self.load:
			print(f"[save] It seems that data was never read or saved. {filename} may contain no data")
			
		stk = np.stack([self.line, self.wave, self.stki, self.stkq, self.stku, self.stkv], axis=1, dtype=np.float64)
		
		np.save(filename, stk)

		del stk  # make the memory free

		return self

	def set_dim(self, num, nw):
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

		self.line = np.zeros(shape=(num, nw), dtype=np.float64)
		self.wave = np.zeros(shape=(num, nw), dtype=np.float64)
		self.stki = np.zeros(shape=(num, nw), dtype=np.float64)
		self.stkq = np.zeros(shape=(num, nw), dtype=np.float64)
		self.stku = np.zeros(shape=(num, nw), dtype=np.float64)
		self.stkv = np.zeros(shape=(num, nw), dtype=np.float64)
		
		self.num = num
		self.nw = nw

		return self


	def write(self, filename, x):
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
			f.write(f" {int(self.line[x,i]):>2} {self.wave[x,i]:>10.4f} {self.stki[x,i]:>14.7E} {self.stkq[x,i]:>14.7E} {self.stku[x,i]:>14.7E} {self.stkv[x,i]:>14.7E}\n")
		f.close()



def read(file1, file2=None):

	pro = Profile(file1, file2)

	return pro

