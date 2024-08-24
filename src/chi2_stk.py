"""
Chi2 Calculation for Stokes Profiles
====================================

This module provides tools to compute, read, and write the chi-squared (χ²) values for Stokes profiles, which are essential in analyzing polarimetric data in astrophysical observations. The χ² statistic is used to measure the goodness of fit between observed and synthetic (modeled) Stokes profiles.

Classes:
--------
- chi2_stk: A class encapsulating the χ² values for Stokes profiles, with methods to compute, read, and write these values.

Functions:
----------
- read_chi2: Reads a binary file containing χ² data and returns a chi2_stk object.
"""

import numpy as np
from scipy.io import FortranFile
import profile_stk as p
import os

def _compute_chi2(obs : np.array, syn : np.array, noise : float, weight : float, num_of_nodes : int, mul=1.0) -> float:
	r"""
	chi2 of one stokes parameter computed with the equation
		$$\chi^2 = \frac{1}{\Lambda} \sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$

	Parameters
	----------
	obs : np.array
		Array with the observations
	syn : np.array
		Array with the inversion model / synthesised Stokes profiles
	noise : float
		Noise of the Stokes parameter
	weight : float
		Weigth of the concerning Stokes parameter
	num_of_nodes : int
		Used number of nodes in the last cycle
	mul : float,optional
		Additional factor to be multiplied to the normalisation

	Returns
	-------
	float
		chi2 value
	"""
	result = 0.0
	
	# Sum over the wavelength points
	for i in range(len(obs)):
		result += (obs[i]-syn[i])*(obs[i]-syn[i])
	
	result = result * weight/noise * weight/noise
	
	result /= (np.float64(mul)*np.float64(len(obs))-np.float64(num_of_nodes))

	return result


def _compute_total_chi2_fov(obs : p.profile_stk, syn : p.profile_stk, n1 : float, n2 : float, n3 : float, n4 : float, w1 : float, w2 : float, w3 : float, w4 : float, num_of_nodes : int) -> float:
	r"""
	$\chi^2$ of all Stokes parameter of the full FOV computed with
		$$\chi^2 = \frac{1}{4n_xn_y\Lambda - F} \sum_x^{n_x}\sum_y^{n_y}\sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$

	Parameters
	----------
	obs : profile_stk
		Observations
	syn : profile_stk
		Synthesis profiles
	n1 : float
		Noise of the parameter stki
	n2 : float
		Noise of the parameter stkq
	n3 : float
 		Noise of the parameter stku
	n4 : float
		Noise of the parameter stkv
	w1 : float
		Weight of the parameter stki
	w2 : float
		Weight of the parameter stkq
	w3 : float
		Weight of the parameter stku
	w4 : float
		Weight of the parameter stkv
 	num_of_nodes : int
 		Number of nodes used for the last inversion "cycle"
		 
    Returns
	-------
	float
	 	reduced chi2 value
	"""
	result = 0.0
	
	for x in range(obs.nx):
		for y in range(obs.ny):
			# Compute chi2 individually for each stokes parameter and add it
			result += _compute_chi2(obs.stki[x][y],syn.stki[x][y],n1,w1,num_of_nodes,4.0*(obs.nx)*(obs.ny))
			result += _compute_chi2(obs.stkq[x][y],syn.stkq[x][y],n2,w2,num_of_nodes,4.0*(obs.nx)*(obs.ny))
			result += _compute_chi2(obs.stku[x][y],syn.stku[x][y],n3,w3,num_of_nodes,4.0*(obs.nx)*(obs.ny))
			result += _compute_chi2(obs.stkv[x][y],syn.stkv[x][y],n4,w4,num_of_nodes,4.0*(obs.nx)*(obs.ny))
		
	
	return result




def _compute_total_chi2(obs : p.profile_stk, syn : p.profile_stk, x : int, y : int, n1 : float, n2 : float, n3 : float, n4 : float, w1 : float, w2 : float, w3 : float, w4 : float, num_of_nodes : int):
	r"""
	$\chi^2$ of all Stokes parameter based on 
		$$\chi^2 = \frac{1}{4\cdot\Lambda - F} \sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$

	Parameters
	----------
	obs : profile_stk
		Observations
	syn : profile_stk
		Synthesis profiles
	x : int
		x position
	y : int
		y position
	n1 : float
		Noise of the parameter stki
	n2 : float
		Noise of the parameter stkq
	n3 : float
 		Noise of the parameter stku
	n4 : float
		Noise of the parameter stkv
	w1 : float
		Weight of the parameter stki
	w2 : float
		Weight of the parameter stkq
	w3 : float
		Weight of the parameter stku
	w4 : float
		Weight of the parameter stkv
 	num_of_nodes : int
 		Number of nodes used for the last inversion "cycle"

    Returns
	-------
	_compute_total_chi2 : float
	 	reduced chi2 value
		 
	"""
	result = 0.0
	
	# Compute chi2 individually for each stokes parameter and add it
	result += _compute_chi2(obs.stki[x][y],syn.stki[x][y],n1,w1,num_of_nodes,4.0)
	result += _compute_chi2(obs.stkq[x][y],syn.stkq[x][y],n2,w2,num_of_nodes,4.0)
	result += _compute_chi2(obs.stku[x][y],syn.stku[x][y],n3,w3,num_of_nodes,4.0)
	result += _compute_chi2(obs.stkv[x][y],syn.stkv[x][y],n4,w4,num_of_nodes,4.0)

	
	return result


class chi2_stk:
	"""
	
	A class to store and manage chi-squared (χ²) values for Stokes profiles. The Stokes profiles
	represent the polarization state of light in astrophysical observations, and the χ² statistic
	is used to assess the goodness of fit between observed and synthetic profiles.

	Parameters
	----------
	total : float
		The total χ² value for the entire map.
	Itot : float
		The total χ² value for Stokes I.
	Qtot : float
		The total χ² value for Stokes Q.
	Utot : float
		The total χ² value for Stokes U.
	Vtot : float
		The total χ² value for Stokes V.
	tot : numpy.ndarray
		A 2D array storing the χ² values for all Stokes parameters combined.
	stki : numpy.ndarray
		A 2D array storing the χ² values for Stokes I.
	stkq : numpy.ndarray
		A 2D array storing the χ² values for Stokes Q.
	stku : numpy.ndarray
		A 2D array storing the χ² values for Stokes U.
	stkv : numpy.ndarray
		A 2D array storing the χ² values for Stokes V.
	nx : int
		The dimension of the data in the x-direction.
	ny : int
		The dimension of the data in the y-direction.
	ns : int
		The number of Stokes parameters (normally 4: I, Q, U, V).
	noise : float
		The noise level for the Stokes parameters.

	Methods
	-------
	compute:
		Computes the χ² values for all pixels across all Stokes parameters and the total χ².
	
	read:
		Reads a binary Fortran file containing χ² data.
	
	write:
		Writes the χ² data to a binary file.
	"""


	def __init__(self, nx : int, ny : int):
		"""
		Initialisation of the class stk_chi2

		Parameters
		----------
		nx : int
			Dimension in x
		ny : int
			Dimension in y

		Returns
		-------
		None
		
		"""
		self.noise = 1e-3 # Weights from SIR
		self.total = 0
		self.Itot = 0
		self.Qtot = 0
		self.Utot = 0
		self.Vtot = 0
		self.tot = np.zeros(shape=(nx,ny))
		self.stki = np.zeros(shape=(nx,ny))
		self.stkq = np.zeros(shape=(nx,ny))
		self.stku = np.zeros(shape=(nx,ny))
		self.stkv = np.zeros(shape=(nx,ny))
		self.nx = nx
		self.ny = ny
		self.ns = 5
		return None
		
	def compute(self, obs : p.profile_stk, syn : p.profile_stk, weights : np.array, num_of_nodes : int) -> None:
		r"""
		Compute the $\chi^2$ of all pixels in all Stokes Parameter and also the total $\chi^2$
		The computation is based on (Borrero et al, 2021) with the equation
		
		$$\chi^2 = \frac{1}{\Lambda - F} \sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$
		
		with $\omega = \frac{w}{\sigma}$ with the used weight $w$ and the noise $\sigma$. In case, there are
		more spectra considered, $\Lambda$ may be multiplied with the number of pixels and number of used Stokes Parameter.
		Obviously, also more summations must be considered.

		Parameters
		----------
		obs : profile_stk
			Observations
		syn : profile_stk
			Synthesis profiles
		weights : np.array
			Weights of the Stokes Parameter
		num_of_nodes : int
			Number of nodes used for the last inversion "cycle"

		Returns
		-------
		None

		Raises
		------
		RuntimeError
			if any dimension of obs or syn is not the same with the class
			
		"""
		if(self.nx != obs.nx):
			raise RuntimeError("[compute] nx of the class and the observations are not the same")
		
		if(self.ny != obs.ny):
			raise RuntimeError("[compute] ny of the class and the observations are not the same")
		
		if(self.ny != syn.ny):
			raise RuntimeError("[compute] nx of the class and the synthesis are not the same")
		
		if(self.ny != syn.ny):
			raise RuntimeError("[compute] ny of the class and the synthesis are not the same")
		
		# Compute chi2 for each pixel
		for x in range(self.nx):
			for y in range(self.ny):
				self.tot[x][y] = _compute_total_chi2(obs, syn, x, y, self.noise, self.noise, self.noise, self.noise,
										 weights[0],weights[1],weights[2], weights[3], num_of_nodes)
				self.stki[x,y] = _compute_chi2(obs.stki[x,y],syn.stki[x,y],self.noise, weights[0], num_of_nodes,1)
				self.stkq[x,y] = _compute_chi2(obs.stkq[x,y],syn.stkq[x,y],self.noise, weights[1], num_of_nodes,1)
				self.stku[x,y] = _compute_chi2(obs.stku[x,y],syn.stku[x,y],self.noise, weights[2], num_of_nodes,1)
				self.stkv[x,y] = _compute_chi2(obs.stkv[x,y],syn.stkv[x,y],self.noise, weights[3], num_of_nodes,1)

		# Compute total chi2
		self.total = _compute_total_chi2(obs, syn, x, y, self.noise, self.noise, self.noise, self.noise, weights[0], weights[1], weights[2], weights[3], num_of_nodes)
		
		self.Itot = np.sum(self.stki) / (np.float64(self.nx) * np.float64(self.ny) * np.float64(obs.nw) - np.float64(num_of_nodes))
		self.Qtot = np.sum(self.stkq) / (np.float64(self.nx) * np.float64(self.ny) * np.float64(obs.nw) - np.float64(num_of_nodes))
		self.Utot = np.sum(self.stku) / (np.float64(self.nx) * np.float64(self.ny) * np.float64(obs.nw) - np.float64(num_of_nodes))
		self.Vtot = np.sum(self.stkv) / (np.float64(self.nx) * np.float64(self.ny) * np.float64(obs.nw) - np.float64(num_of_nodes))

		return



	def read(self, fname : str, fmt_type=np.float32):
		"""
		Read a binary fortran file

		Parameters
		----------
		fname : str
			Name of the binary file
		fmt_type : type
			Type of the fortran file
		
		Returns
		-------
		chi2_stk
			Class with the read data

		"""
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

		if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-110904.) > 1.e-3) ):
			print('Something is wrong with the file %s' % fname)
			print('\tIs it a 2D chi2 map file?')
			return np.nan

		if (np.abs(vers-3.) < 1.e-3):
			#VERSION 3
			dims = first_rec[5:]
			nx, ny, _, ns = np.int16(dims)
			self.nx = nx
			self.ny = ny
			self.ns = ns
		
			# Read total chi2 for the full map
			tmp = f.read_record(dtype=fmt_type)
			array = np.zeros(5, dtype=np.float32)
			array[0:tmp.size] = tmp*1.
			
			self.total = array[0] # Total chi2 value
			self.Itot = array[1] # Total chi2 value
			self.Qtot = array[2] # Total chi2 value
			self.Utot = array[3] # Total chi2 value
			self.Vtot = array[4] # Total chi2 value

			# Read the data for chi2 total, stki, stkq, stku and stkv
			array = np.zeros(np.int64((ns*1.)*nx*ny), dtype=np.float32)
		
			tmp = f.read_record(dtype=fmt_type)
			array[0:tmp.size] = tmp*1.

			array = np.moveaxis(array.reshape(nx,ny,5)\
        				, [0,1,2], [1,2,0]) * 1.

			self.tot = array[0,:,:]
			self.stki = array[1,:,:]
			self.stkq = array[2,:,:]
			self.stku = array[3,:,:]
			self.stkv = array[4,:,:]

			del array

		else:
			print('Version %i of chi2 file is not supported' % np.int(vers))
			return np.nan

		f.close()

		return self
	
	def write(self, fname : str, fmt_type=np.float32) -> None:
		"""
		Write the data to a binary file

		Parameters
		----------
		fname : str
			File name
		fmt_type : type
			Type of the fortran file

		"""
		f = FortranFile(fname, 'w')
		metadata = np.zeros(9,dtype=np.float32)

		# Fill Metadata
		SS = 3000
		VV = 3

		metadata[0] = (SS + VV)
		metadata[1] = (SS - VV)
		metadata[2] = (110904.0)
		metadata[3] = 2.0
		metadata[4] = 3.0
		metadata[5] = self.nx
		metadata[6] = self.ny
		metadata[7] = 0.0
		metadata[8] = self.ns		

		f.write_record(metadata)

		
		# Write total chi2 for the full map
		total = np.zeros(5,dtype=np.float32)
		total[0] = self.total
		total[1] = self.Itot
		total[2] = self.Qtot
		total[3] = self.Utot
		total[4] = self.Vtot

		f.write_record(total)
		

		# Write the data for chi2 total, stki, stkq, stku and stkv
		array = np.zeros(np.int64((self.ns*1.)*self.nx*self.ny), dtype=np.float32)
		index = 0
		for j in range(self.nx):
			for i in range(self.ny):
				array[index] = self.tot[j][i]
				array[index+1] = self.stki[j][i]
				array[index+2] = self.stkq[j][i]
				array[index+3] = self.stku[j][i]
				array[index+4] = self.stkv[j][i]
				index += 5
		
		f.write_record(array)
	
		f.close()

		return self
	
def read_chi2(filename : str, fmt_type = np.float32):
	"""
	Reads a chi2 file and loads the data

	Parameters
	----------
	filename : str
		Filename of the binary file
	fmt_type : type
		Type of the binary file (only np.float32 implemented)

	Returns
	-------
	read_chi2 : chi2_stk
		Class with the read binary file
		
	Raises
	------
	FileExistsError
		if first file does not exist
	"""
	if not os.path.exists(filename):
		raise FileExistsError("[read_chi2] " + filename + " does not exist.")
	
	f = FortranFile(filename, 'r')
	first_rec = f.read_record(dtype=fmt_type)

	posv = first_rec[0]
	negv = first_rec[1]
	fid = first_rec[2]

	vers = (posv - negv) / 2.

	if (np.abs(fid-110904.) > 1.e-3):
		print('Something is wrong with the file %s' % filename)
		print('\tIs it a 2D chi2 map file?')
		return np.nan

	if (np.abs(vers-3.) < 1.e-3):
		#VERSION 3
		dims = first_rec[5:]
		nx, ny, _, ns = np.int16(dims)
		chi2 = chi2_stk(nx,ny)
		chi2.read(filename, fmt_type)
		
		return chi2
	else:
		print('Version %i of chi2 file is not supported' % np.int(vers))
		return np.nan
	
