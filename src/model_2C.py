import numpy as np 
import os


class Model:
	"""
	Class containing the models of a simulation

	Variables are:
		- log_tau
		- T
		- Pe
		- vmicro
		- B
		- vlos
		- gamma
		- phi
		- z
		- pg
		- rho


	There are several functions:
		- read: Read a numpy file containing models and stores it into the class variables
		- read_results: Reads the results from the inversion of SIR
		- write: Writes a Model to a SIR readable file
		- save: Saves the models as a numpy file
		- set_limit: Cuts the data to a specific log_tau value (used for plotting)


	"""
	def __init__(self, filename=None, nx = 0, ny = 0, npar = 0, filename_fill=None):
		"""
		Initialisation of the class with the models
		
		Parameter
		---------
		filename : str, optional
			File name of the model to be loaded. Default: None (= no reading)
		filename_fill : str, optional
			File name of the filling factor
		
		"""
		self.full = False		# Determines whether also z, pg, and rho exists
		self.load = False		# Determines whether data is already loaded or not
		self.num = 0			# Number of Models saved
		self.nx = nx				# Points in x
		self.ny = ny				# Points in y
		self.npar = npar			# Points in log tau
		if filename is not None:
			mods = np.load(filename).astype(np.float64)
			self.log_tau = mods[:, :, 0, :]
			self.T = mods[:, :, 1, :]
			self.Pe = mods[:, :, 2, :]
			self.vmicro = mods[:, 3, :]
			self.B = mods[:, :, 4, :]
			self.vlos = mods[:, :, 5, :]/1e5
			self.gamma = mods[:, :, 6, :]
			self.phi = mods[:, :, 7, :]
			if len(mods) > 8:
				self.z = mods[:, :, 8, :]
				self.Pg = mods[:, :, 9, :]
				self.rho = mods[:, :, 10, :]
				self.full = True
			else:
				self.z = np.zeros(shape=self.B.shape, dtype=np.float64)
				self.Pg = np.zeros(shape=self.B.shape, dtype=np.float64)
				self.rho = np.zeros(shape=self.B.shape, dtype=np.float64)
				self.full = False
			self.load = True
			self.nx = self.log_tau.shape[0]
			self.ny = self.log_tau.shape[1]
			self.npar = self.log_tau.shape[2]
			if filename_fill is not None:
				self.fill = np.load(filename_fill).astype(np.float64)
			else:
				self.fill = np.zeros(shape=(nx,ny), dtype=np.float64)
		else:
			self.log_tau = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.T = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.Pe = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.vmicro = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.B = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.vlos = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.gamma = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.phi = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.z = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.Pg = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.rho = np.zeros(shape=(nx,ny,npar), dtype=np.float64)
			self.fill = np.zeros(shape=(nx,ny), dtype=np.float64)

	def correct_phi(self):
		"""
		SIR can give you any value for the azimuth, e.g. -250 deg, and therefore, to compute the standard deviation
		those values should be corrected so that I have values from 0 to 360 degrees.

		Parameter
		---------
		None
		
		Return
		------
		class
		"""
		# No correction needed
		if len(self.phi[self.phi < 0]) == len(self.phi[self.phi > 180]) == 0:
			return self

		# Correct for values to be between 0 and 180 (we cannot measure differences between 45 and 225 deg)
		# Note that negative values are also correctly transformed (e.g. -175 % 180 = 5)
		self.phi = self.phi % 180

	def get_attribute(self, string):
		"""
		Returns a specific physical parameter. This can be used if ones only wants a specific
		parameter depending on a string/input

		Parameter
		---------
		string : str
			Determines which physical parameter is returned
			Options are: tau, T, Pe, vmicro, B, vlos, gamma, phi, z, Pg, rho, nx, ny, npar, fill

		Return
		------
		The wished physical parameter
		"""
		
		if string == "tau":
			return self.log_tau
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

	def read(self, filename, filename_fill = ""):
		"""
		Reads a numpy file and stores the parameters
		
		Parameter
		---------
		filename : str
			File name of the model to be loaded
		filename_fill : str
			Filename with the filling factors
		
		"""
		mods = np.load(filename).astype(np.float64)
		self.log_tau = mods[:, :, 0, :]
		self.T = mods[:, :, 1, :]
		self.Pe = mods[:, :, 2, :]
		self.vmicro = mods[:, 3, :]
		self.B = mods[:, :, 4, :]
		self.vlos = mods[:, :, 5, :] / 1e5
		self.gamma = mods[:, :, 6, :]
		self.phi = mods[:, :, 7, :]
		if len(mods) > 8:
			self.z = mods[:, :, 8, :]
			self.Pg = mods[:, :, 9, :]
			self.rho = mods[:, :, 10, :]
			self.full = True
		else:
			self.z = np.zeros(shape=self.B.shape, dtype=np.float64)
			self.Pg = np.zeros(shape=self.B.shape, dtype=np.float64)
			self.rho = np.zeros(shape=self.B.shape, dtype=np.float64)
			self.full = False
		self.load = True
		self.nx = self.log_tau.shape[0]
		self.ny = self.log_tau.shape[1]
		self.npar = self.log_tau.shape[2]

		if filename_fill != "":
			self.fill = np.load(filename).astype(np.float64)
		else:
			self.fill = np.zeros(shape=(self.nx,self.ny), dtype=np.float64)
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
		temp = np.loadtxt(f"{os.path.join(path,task['folders'][0])}/{filename}", skiprows=1).transpose()
		self.npar = len(temp[0])

		self.log_tau = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.T = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.Pe = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.vmicro = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.B = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.vlos = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.gamma = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.phi = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.z = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.Pg = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.rho = np.zeros(shape=(nx, ny, self.npar), dtype=np.float64)
		self.fill = np.zeros(shape=(nx, ny), dtype=np.float64)
		
		del temp  # Remove the temporary array from the memory

		for i in range(len(task['x'])):
			file = f"{os.path.join(path,task['folders'][i])}/{filename}"
			data = np.loadtxt(file, skiprows=1, dtype=np.float64).transpose()
			header = np.genfromtxt(file,max_rows=1)
			x, y = task['x'][i]-task['x'][0], task['y'][i]-task['y'][0]
			self.log_tau[x, y] = data[0]
			self.T[x, y] = data[1]
			self.Pe[x, y] = data[2]
			self.vmicro[x, y] = data[3]
			self.B[x, y] = data[4]
			self.vlos[x, y] = data[5]/1e5
			self.gamma[x, y] = data[6]
			self.phi[x, y] = data[7]
			if len(data) > 8:
				self.z[x, y] = data[8]
				self.Pg[x, y] = data[9]
				self.rho[x, y] = data[10]
			self.fill[x,y] = float(header[1])
		self.load = True
		self.full = True
		
		return self

	def save(self, filename, fill=False):
		"""
		Save the data in the class as a numpy file

		Parameters
		----------
		filename : str
			Filename of the numpy file
		fill : bool, optional
			Also save filling factor. Default: False

		"""
		if not self.load:
			print(f"[save] It seems that data was never read or saved. {filename} max contain no data")

		# Consistency check for vlos
		if np.mean(self.vlos) > 1e4:
			print("[save] It seems that the data in the class are in cm/s and not in km/s as expected. It is corrected.")
			print("       If this is not wished or wrong, look at the function save in model.py.")
			self.vlos = self.vlos/1e5
		
		if (self.rho == np.zeros(shape=self.log_tau.shape)).all(): # No information in z, rho and Pg and thereofre do not save it
			models = np.stack([self.log_tau, self.T, self.Pe, self.vmicro, self.B, self.vlos*1e5, self.gamma,
									self.phi], axis=2, dtype=np.float64)
		else:
			models = np.stack([self.log_tau, self.T, self.Pe, self.vmicro, self.B, self.vlos*1e5, self.gamma,
									self.phi, self.z, self.Pg, self.rho], axis=2, dtype=np.float64)
		np.save(filename, models)
		if fill:
			np.save(filename.replace(".npy", "_fill.npy"), self.fill)

		del models  # make the memory free

		return self

	def set_dim(self, nx, ny, npars):
		"""
		Sets the dimensions if no data is loaded yet

		Parameter
		---------
		nx : int
			Number of Models in x
		ny : int
			Number of Models in y
		npars : int
			Number of values per physical parameter

		"""
		if self.load:
			print("[set_dim] Data is already loaded and therefore dimensions cannot be set.")
			return self

		self.log_tau = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.T = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.Pe = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.vmicro = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.B = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.vlos = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.gamma = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.phi = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.z = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.Pg = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.rho = np.zeros(shape=(nx, ny, npars), dtype=np.float64)
		self.fill = np.zeros(shape=(nx, ny), dtype=np.float64)
		self.nx = nx * 1
		self.nx = ny * 1

		return self

	def set_limit(self, lim):
		"""
		Cuts the arrays to a specific log tau value (should only be used for plotting)

		Parameter
		---------
		lim : int
			log tau value to which the values are cut

		"""
		if self.log_tau[0, 0, -1] > lim:
			ind = np.argmin(np.abs(self.log_tau[0, 0, :] - lim)) + 1
			self.log_tau = self.log_tau[:, :, 0:ind]
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

	def write(self, filename, header, x, y):
		"""
		Write a model with the given data in a specific format.

		Parameter
		---------
		filename : string
			Name of the saved file
		Header : string
			Header of the model
		x : int
			Integer determining which model
		y : int
			Integer determining which model

		Return
		------
		None
		"""
		if x > self.nx:
			print(f"[write] Error in writing the model in model.py. x = {x} is bigger than the dimension {self.nx} of the array")
		elif x > self.nx:
			print(f"[write] Error in writing the model in model.py. y = {y} is bigger than the dimension {self.ny} of the array")
		if x > self.nx or y > self.ny:
			return self
		
		# Correct that vlos is most likely wrong and in cm/s instead of km/s
		if np.max(self.vlos) > 1.0e4:
			self.vlos = self.vlos / 1.0e5 

		f = open(filename, 'w')
		f.write(f"{header}\n")
		if self.full:
			for n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11 in zip(self.log_tau[x, y], self.T[x, y], self.Pe[x, y],
																	self.vmicro[x, y], self.B[x, y], self.vlos[x, y]*1.0e5,
																	self.gamma[x, y], self.phi[x, y], self.z[x, y],
																	self.Pg[x, y], self.rho[x, y]):
				f.write(f" {n1:>7.4f} {n2:>7.1f} {n3:>12.5E} {n4:>10.3E} {n5:>11.4E} {n6:>11.4E} {n7:>11.4E} {n8:>11.4E} {n9:>11.4E} {n10:>11.4E} {n11:>11.4E}\n")

		else:
			for n1, n2, n3, n4, n5, n6, n7, n8 in zip(self.log_tau[x, y], self.T[x, y], self.Pe[x, y], self.vmicro[x, y],
														self.B[x, y], self.vlos[x, y] * 1.0e5, self.gamma[x, y],
														self.phi[x, y]
													):
				f.write(f" {n1:>7.4f} {n2:>7.1f} {n3:>12.5E} {n4:>10.3E} {n5:>11.4E} {n6:>11.4E} {n7:>11.4E} {n8:>11.4E}\n")

		return self

	
class Error(Model):
	"""
	Class containing the errors of a simulation

	Variables are:
		- log_tau
		- T
		- Pe
		- vmicro
		- B
		- vlos
		- gamma
		- phi
		- z
		- pg
		- rho


	This class contains the same functions as the class model but adapted to the error

	"""
	pass

	def save(self, filename):
		"""
		Save the data in the class as a numpy file

		Parameters
		----------
		filename : str
			Filename of the numpy file
		
		"""
		if not self.load:
			print(f"[save] It seems that data was never read or saved. {filename} max contain no data")
			
		if (self.rho == np.zeros(shape=self.log_tau.shape)).all(): # No information in z, rho and Pg and thereofre do not save it
			models = np.stack([self.log_tau, self.T, self.Pe, self.vmicro, self.B, self.vlos*1e5, self.gamma,
									self.phi], axis=2, dtype=np.float64)
		else:
			models = np.stack([self.log_tau, self.T, self.Pe, self.vmicro, self.B, self.vlos*1e5, self.gamma,
									self.phi, self.z, self.Pg, self.rho], axis=2, dtype=np.float64)
		np.save(filename, models)
		

		del models  # make the memory free

		return self
