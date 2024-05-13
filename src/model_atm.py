"""
Model
=====

Class Model with all the tools to read and write the models. The class contains all the information from the model.

"""
import numpy as np 
import os
from scipy.io import FortranFile

class Model:
	"""
	Class containing the models of a simulation

	Parameters
	----------
	tau : numpy.array
		Array with the log tau
	T : numpy.ndarray
		Temperature in K
	Pe : numpy.ndarray
		Electron Pressure in dyn/cm$^2$
	vmicro : numpy.ndarray
		Microturbulence Velocity in cm/s
	B : numpy.ndarray
		Magnetic Field Strength in G
	vlos : numpy.ndarray
		Line-of-Sight Velocity in km/s
	gamma : numpy.ndarray
		Inclination in deg
	phi : numpy.ndarray
		Azimuth in deg
	z : numpy.ndarray
		z in km
	Pg : numpy.ndarray
		Gas pressure in dyn/cm$^2$
	rho : numpy.ndarray
		Density in g/cm$^3$
	vmacro : numpy.array
		Macroturbulence Velocity in km/s
	fill : numpy.array
		Filling factor
	stray_light : numpy.array
		Stray Light Factor in Percent

	Returns
	-------
	None


	"""
	def __init__(self, nx = 0, ny = 0, nval = 0):
		"""

		Initialisation of the class with the models.
		
		Parameters
		----------
		nx : int, optional
			Dimension in x. Default: 0.
		ny : int, optional
			Dimension in y. Default: 0.
		nval : int, optional
			Dimension in "z". Default: 0.
		
		Returns
		-------
		None
		
		"""
		self.full = False			# Determines whether also z, pg, and rho exists
		self.load = False			# Determines whether data is already loaded or not
		self.nx = nx				# Points in x
		self.ny = ny				# Points in y
		self.nval = nval			# Points in log tau
		self.tau = np.zeros(shape=(nval), dtype=np.float64)
		self.T = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.Pe = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.vmicro = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.B = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.vlos = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.gamma = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.phi = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.z = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.Pg = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.rho = np.zeros(shape=(nx,ny,nval), dtype=np.float64)
		self.vmacro = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.fill = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.stray_light = np.zeros(shape=(nx, ny), dtype=np.float64)

	def correct_phi(self):
		"""

		SIR can give you any value for the azimuth, e.g. -250 deg, and therefore, to compute the standard deviation
		those values should be corrected so that I have values from 0 to 360 degrees.


		Parameters
		----------
		None
		
		Returns
		-------
		None

		"""
		# No correction needed
		if len(self.phi[self.phi < 0]) == len(self.phi[self.phi > 180]) == 0:
			return self

		# Correct for values to be between 0 and 180 (we cannot measure differences between 45 and 225 deg)
		# Note that negative values are also correctly transformed (e.g. -175 % 180 = 5)
		self.phi = self.phi % 180

	def cut_to_map(self, Map):
		"""
		Cut the data to a map [xmin, xmax, ymin, ymax]

		Parameters
		----------
		Map : list
			List with the ranges in pixel in x and y direction

		Returns
		-------
		None

		"""
		
		self.nx = Map[1]-Map[0]+1
		self.ny = Map[3]-Map[2]+1

		self.T = self.T[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.Pg = self.Pg[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.vmicro = self.vmicro[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.B = self.B[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.vlos = self.vlos[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.gamma = self.gamma[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.phi = self.phi[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.z = self.z[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.rho = self.rho[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.Pe = self.Pe[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.vmacro = self.vmacro[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.fill = self.fill[Map[0]:Map[1]+1, Map[2]:Map[3]+1]
		self.stray_light = self.stray_light[Map[0]:Map[1]+1, Map[2]:Map[3]+1]

		
		
		return self

	def get_attribute(self, string):
		"""
		Returns a specific physical parameter. This can be used if ones only wants a specific
		parameter depending on a string/input.

		Parameters
		----------
		string : str
			Determines which physical parameter is returned
			Options are: tau, T, Pe, vmicro, B, vlos, gamma, phi, z, Pg, rho, nx, ny, npar, fill

		Returns
		-------
		out : float
			The wished physical parameter

		"""
		
		if string == "tau":
			return self.tau
		if string == "T":
			return self.T
		if string == "Pe":
			return self.Pe
		if string == "vmicro":
			return self.vmicro
		if string == "B":
			return self.B
		if string == "vlos":
			return self.vlos
		if string == "gamma":
			return self.gamma
		if string == "phi":
			return self.phi
		if string == "z":
			return self.z
		if string == "Pg":
			return self.Pg
		if string == "rho":
			return self.rho
		if string == "nx":  # Number of Models in x
			return self.nx
		if string == "ny":  # Number of Models in y
			return self.ny
		if string == "npar":
			return self.npar
		if string == "fill":
			return self.fill
		return None

	def interp(self, new_tau):
		"""
		Interpolate the model to a new log tau scale.

		Parameters
		----------
		new_tau : numpy array
			New log tau scale
		
		Returns
		-------
		None

		"""
		# Perform some checks
		if self.nval != len(new_tau):
			print(f"[interp] Error as the shape is incorecct ({len(new_tau)} vs. {self.nval})")
			return self
		if new_tau[-1] < self.tau[-1]:
			print(f"[interp] Error as the scale exceeds the borders of the model ({new_tau[-1]} vs. {self.tau[-1]})")
			return self
		if new_tau[0] > self.tau[0]:
			print(f"[interp] Error as the scale exceeds the borders of the model ({new_tau[0]} vs. {self.tau[0]})")
			return self
		
		for x in self.nx:
			for y in self.y:
				self.T[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.T[x,y]))
				self.Pe[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.Pe[x,y]))
				self.vmicro[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.vmicro[x,y]))
				self.B[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.B[x,y]))
				self.vlos[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.vlos[x,y]))
				self.gamma[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.gamma[x,y]))
				self.phi[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.phi[x,y]))
				if self.full or np.any(self.z):
					self.z[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.z[x,y]))
					self.Pg[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.Pg[x,y]))
					self.rho[x,y] = np.interp(new_tau, np.flip(self.tau), np.flip(self.rho[x,y]))
		self.tau = new_tau

		return self

	def read(self, fname, fmt_type=np.float64):
		"""
		Reads a binary model file.

		Parameters
		----------
		fname : str
			Filename
		fmt_type : type, optional
			Type used to laod data. Default: np.float64

		Returns
		-------
		None

		"""
		f = FortranFile(fname, 'r')
		first_rec = f.read_record(dtype=np.float32) # Must be float32 for compatibility reasons

		posv = first_rec[0]
		negv = first_rec[1]
		fid = first_rec[2]
		nrec = first_rec[3]
		ndims = first_rec[4]

		medv = (posv + negv) / 2.
		verp = posv - medv
		vern = negv - medv
		vers = (posv - negv) / 2.

		if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-230904.) > 1.e-3) ):
			print('Something is wrong with the file %s' % fname)
			print('\tIs it a 3D SIR MIG model file?')
			return np.nan

		if (np.abs(vers-3.) < 1.e-3):
			#VERSION 3
			dims = first_rec[5:]
			npar, nx, ny, nval = np.int16(dims)
			self.nx = nx
			self.ny = ny
			self.npar = npar
			self.nval = nval

			# Read log tau values
			tmp = f.read_record(dtype=fmt_type)
			self.tau = tmp * 1.
			
			# Read vmacro
			tmp = f.read_record(dtype=fmt_type)
			array = tmp * 1.
			self.vmacro = array.reshape(self.nx, self.ny)

			# Read filling factor
			tmp = f.read_record(dtype=fmt_type)
			array = tmp * 1.
			self.fill = array.reshape(self.nx, self.ny)

			# Read fstray light factor
			tmp = f.read_record(dtype=fmt_type)
			array = tmp * 1.
			self.stray_light = array.reshape(self.nx, self.ny)


			# Read model parameter
			tmp = f.read_record(dtype=fmt_type)
			array = tmp.reshape(self.nx,self.ny,self.npar, self.nval) * 1.
			self.T = array[:,:,0]
			self.Pe = array[:,:,1]
			self.vmicro = array[:,:,2]
			self.B = array[:,:,3]
			self.vlos = array[:,:,4]/1e5
			self.gamma = array[:,:,5]
			self.phi = array[:,:,6]
			if npar > 9:
				self.z = array[:,:,7]
				self.Pg = array[:,:,8]
				self.rho = array[:,:,9]
				self.full = True
			else:
				self.z = np.zeros(shape=self.T.shape)
				self.Pg = np.zeros(shape=self.T.shape)
				self.rho = np.zeros(shape=self.T.shape)
				self.full = False

			self.load = True

			del array

		else:
			print('Version %i of model file is not supported' % np.int(vers))
			return np.nan

		f.close()

		return self
	
	def read_mod(self, fname, fmt_type=np.float64):
		"""
		Reads a model file from SIR and assigns the value to the class.
		This class will then create the class with nx = ny = 1. 

		Parameters
		----------
		fname : str
			Name of the .mod file

		Returns
		-------
		None

		"""
		if self.load:
			print("[read_mod] Data is already loaded! Change load of this class to False.")
			return self
		
		# Load Data
		data = np.loadtxt(fname, skiprows=1).transpose()
		header = np.genfromtxt(fname,max_rows=1)

		# Set Dimensions
		self.set_dim(1,1,len(data[0]))

		# Store data
		self.tau	= data[0]
		self.T[0,0]		= data[1]
		self.Pe[0,0]	 	= data[2]
		self.vmicro[0,0] 	= data[3]
		self.B[0,0]	 	= data[4]
		self.vlos[0,0]     	= data[5]
		self.gamma[0,0]  	= data[6]
		self.phi[0,0]	= data[7]
		if len(data) > 8:
			self.z[0,0]   = data[8]
			self.Pg[0,0]  = data[9]
			self.rho[0,0] = data[10]
			self.full = True
		else:
			self.full = False

		self.load = True
		
		self.vmacro[0,0] = header[0]
		self.fill[0,0] = header[1]
		self.stray_light[0,0] = header[2]

		del data


		return self


	def read_results(self, task, filename, path, nx, ny):
		"""
		Reads all the errors from the inversion
		
		Parameters
		----------
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
		

		Returns
		-------
		None

		"""	
		# Set the dimensions
		self.nx = nx
		self.ny = ny
		temp = np.loadtxt(f"{os.path.join(path,task['folders'][0])}/{filename}", skiprows=1).transpose()
		self.nval = len(temp[0])

		self.tau = np.zeros(shape=(self.nval), dtype=np.float64)
		self.T = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.Pe = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.vmicro = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.B = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.vlos = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.gamma = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.phi = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.z = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.Pg = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.rho = np.zeros(shape=(nx, ny, self.nval), dtype=np.float64)
		self.vmacro = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.fill = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.stray_light = np.zeros(shape=(nx, ny), dtype=np.float64)
		
		del temp  # Remove the temporary array from the memory

		for i in range(len(task['x'])):
			file = f"{os.path.join(path,task['folders'][i])}/{filename}"
			data = np.loadtxt(file, skiprows=1, dtype=np.float64).transpose()
			header = np.genfromtxt(file,max_rows=1)
			x, y = task['x'][i]-task['x'][0], task['y'][i]-task['y'][0]
			if i == 0:
				self.tau = data[0]
			self.T[x, y] = data[1]
			self.Pe[x, y] = data[2]
			self.vmicro[x, y] = data[3]
			self.B[x, y] = data[4]
			self.vlos[x, y] = data[5]/1e5
			self.gamma[x, y] = data[6]
			self.phi[x, y] = data[7]
			if len(data) > 9:
				self.z[x, y] = data[8]
				self.Pg[x, y] = data[9]
				self.rho[x, y] = data[10]
				self.full = False
			else:
				self.full = True
			self.vmacro[x,y] = float(header[0])
			self.fill[x,y] = float(header[1])
			self.stray_light[x,y] = float(header[2])

		self.load = True
		
		
		return self

	
	def set_dim(self, nx, ny, nval):
		"""
		Sets the dimensions if no data is loaded yet

		Parameters
		----------
		nx : int
			Number of Models in x
		ny : int
			Number of Models in y
		nval : int
			Number of values per physical parameter

		Returns
		-------
		None
		"""
		if self.load:
			print("[set_dim] Data is already loaded and therefore dimensions cannot be set.")
			return self

		self.tau = np.zeros(shape=(nval), dtype=np.float64)
		self.T = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.Pe = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.vmicro = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.B = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.vlos = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.gamma = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.phi = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.z = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.Pg = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.rho = np.zeros(shape=(nx, ny, nval), dtype=np.float64)
		self.vmacro = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.fill = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.stray_light = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.nx = nx * 1
		self.nx = ny * 1
		self.nval = nval * 1

		return self

	def set_limit(self, lim):
		"""

		Cuts the arrays to a specific log tau value (should only be used for plotting).

		Parameters
		----------
		lim : int
			log tau value to which the values are cut

		Returns
		-------
		None

		"""
		if self.tau[-1] > lim:
			ind = np.argmin(np.abs(self.tau - lim)) + 1
			self.tau = self.tau[0:ind]
			self.T = self.T[:, :, 0:ind]
			self.Pe = self.Pe[:, :, 0:ind]
			self.vmicro = self.vmicro[:, :, 0:ind]
			self.B = self.B[:, :, 0:ind]
			self.vlos = self.vlos[:, :, 0:ind]
			self.gamma = self.gamma[:, :, 0:ind]
			self.phi = self.phi[:, :, 0:ind]
			self.z = self.z[:, :, 0:ind]
			self.Pg = self.Pg[:, :, 0:ind]
			self.rho = self.rho[:, :, 0:ind]

		return self

	def write(self, fname, fmt_type=np.float64):
		"""
		Write data into a binary fortran file

		Parameters
		----------
		fname : str
			File name 
		fmt_type : type
			type which is used to save it => numpy.float64 used
		
		Returns
		-------
		None

		"""
		if (fmt_type != np.float64):
			print('Not implemented yet!')
			return np.nan

		if not self.load:
			print(f"[write] It seems that data was never read or saved. {fname} may contain no data")

		# Consistency check for vlos
		if np.mean(np.abs(self.vlos)) > 50:
			self.vlos = self.vlos/1e5

		f = FortranFile(fname, 'w')
		
		# FIRST, WE WRITE A FIRST RECORD WITH THE MANDATORY DATA:
		towrite = np.zeros(9, dtype=fmt_type)
		
		v = 3

		towrite[0] = 3000 + v
		towrite[1] = 3000 - v
		towrite[2] = 230904. * 1.
		towrite[3] = 1. #NREC
		towrite[4] = 4.
		if not np.any(self.rho):
			towrite[5] = 7 * 1.
		else:
			towrite[5] = 10 * 1.
		towrite[6] = self.nx * 1.
		towrite[7] = self.ny * 1.
		towrite[8] = self.T.shape[2] * 1. # Number of values

		f.write_record(np.float32(towrite))
		
		f.write_record(np.float64(self.tau))

		f.write_record(np.float64(self.vmacro.flatten()))
		f.write_record(np.float64(self.fill.flatten()))
		f.write_record(np.float64(self.stray_light.flatten()))

		if np.any(self.rho):
			array = np.zeros(shape=(self.nx,self.ny,10,self.nval))
		else:
			array = np.zeros(shape=(self.nx,self.ny,7,self.nval))


		array[:,:,0] = self.T
		array[:,:,1] = self.Pe
		array[:,:,2] = self.vmicro
		array[:,:,3] = self.B
		array[:,:,4] = self.vlos*1e5
		array[:,:,5] = self.gamma
		array[:,:,6] = self.phi
		if np.any(self.rho):
			array[:,:,7] = self.z
			array[:,:,8] = self.Pg
			array[:,:,9] = self.rho

		f.write_record(np.float64(array.flatten()))

		del array

		f.close()

		return self

	def write_model(self, filename, x, y):
		"""
		Write a model with the given data in a specific format.

		Parameters
		----------
		filename : str
			Name of the saved file
		x : int
			Integer determining which model
		y : int
			Integer determining which model

		Returns
		-------
		None

		"""
		if x > self.nx:
			print(f"[write] Error in writing the model in model.py. x = {x} is bigger than the dimension {self.nx} of the array")
		elif x > self.nx:
			print(f"[write] Error in writing the model in model.py. y = {y} is bigger than the dimension {self.ny} of the array")
		if x > self.nx or y > self.ny:
			return self
		
		# Correct that vlos is most likely wrong and in cm/s instead of km/s
		if np.mean(np.abs(self.vlos)) > 50:
			self.vlos = self.vlos / 1.0e5 

		f = open(filename, 'w')
		f.write(f"   {self.vmacro[x,y]}  {self.fill[x,y]}  {self.stray_light[x,y]}\n")
		if self.full or np.any(self.rho):
			for n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11 in zip(self.tau, self.T[x, y], self.Pe[x, y],
																	self.vmicro[x, y], self.B[x, y], self.vlos[x, y]*1.0e5,
																	self.gamma[x, y], self.phi[x, y], self.z[x, y],
																	self.Pg[x, y], self.rho[x, y]):
				f.write(f" {n1:>7.4f} {n2:>7.1f} {n3:>12.5E} {n4:>10.3E} {n5:>11.4E} {n6:>11.4E} {n7:>11.4E} {n8:>11.4E} {n9:>11.4E} {n10:>11.4E} {n11:>11.4E}\n")

		else:
			for n1, n2, n3, n4, n5, n6, n7, n8 in zip(self.tau, self.T[x, y], self.Pe[x, y], self.vmicro[x, y],
														self.B[x, y], self.vlos[x, y] * 1.0e5, self.gamma[x, y],
														self.phi[x, y]
													):
				f.write(f" {n1:>7.4f} {n2:>7.1f} {n3:>12.5E} {n4:>10.3E} {n5:>11.4E} {n6:>11.4E} {n7:>11.4E} {n8:>11.4E}\n")

		return self

def read_model(filename):
	mod = Model(0,0,0)
	if ".mod" in filename:
		mod.read_mod(filename)
	else:
		mod.read(filename)

	return mod
	
