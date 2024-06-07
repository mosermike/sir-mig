"""

Chi2
====

Class chi2_stk with all the tools to compute, read and write chi2 value of Stokes Profiles.

"""
import numpy as np
from scipy.io import FortranFile

def _compute_chi2(obs, syn, noise, weight, num_of_nodes, mul=1.):
	r"""
	chi2 of one stokes parameter computed with the equation
		$$\chi^2 = \frac{1}{\Lambda} \sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$


	"""
	result = 0.0
	
	# Sum over the wavelength points
	for i in range(len(obs)):
		result += (obs[i]-syn[i])*(obs[i]-syn[i])
	
	result = result * weight/noise * weight/noise
	
	result /= (mul*len(obs)-num_of_nodes)

	return result


def _compute_total_chi2_fov(obs, syn, n1, n2, n3, n4, w1, w2, w3, w4, num_of_nodes):
	r"""
	$\chi^2$ of all Stokes parameter of the full FOV computed with
		$$\chi^2 = \frac{1}{4n_xn_y\Lambda - F} \sum_x^{n_x}\sum_y^{n_y}\sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$

	Parameter
	---------
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
		 
    Return
	------
	out : float
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




def _compute_total_chi2(obs, syn, x, y, n1, n2, n3, n4,	w1, w2, w3, w4, num_of_nodes):
	r"""
	$\chi^2$ of all Stokes parameter based on 
		$$\chi^2 = \frac{1}{4\cdot\Lambda - F} \sum_i^\Lambda (I^{obs}_i - I^{syn}_i)^2 \cdot \omega^2$$

	Parameter
	---------
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

    Return
	------
	out : float
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
	r"""
	Class containing the $\chi^2$-values

	Parameters
	----------
	total : float
		Total chi2 of the full map
	Itot : float
		Total chi2 for Stokes I
	Qtot : float
		Total chi2 for Stokes Q
	Utot : float
		Total chi2 for Stokes U
	Vtot : float
		Total chi2 for Stokes V
	tot : numpy.ndarray
		chi2 for all 4 Stokes Parameter
	stki : numpy.ndarray
		chi2 for Stokes Parameter I
	stkq : numpy.ndarray
		chi2 for Stokes Parameter Q
	stku : numpy.ndarray
		chi2 for Stokes Parameter U
	stkv : numpy.ndarray
		chi2 for Stokes Parameter V
	nx : int
		Dimension in x
	ny : int
		Dimension in y
	ns : int
		Number of Stokes Parameter
	noise : float
		Noise level of the Stokes Parameter

	"""
	def __init__(self, nx, ny):
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
		
	def compute(self, obs, syn, weights, num_of_nodes):
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
		weights : array
			Weights of the Stokes Parameter
		num_of_nodes : int
			Number of nodes used for the last inversion "cycle"

			
		"""
		self.nx = obs.nx
		self.ny = obs.ny
		self.tot = np.zeros(shape=(self.nx,self.ny))
		self.stki = np.zeros(shape=(self.nx,self.ny))
		self.stkq = np.zeros(shape=(self.nx,self.ny))
		self.stku = np.zeros(shape=(self.nx,self.ny))
		self.stkv = np.zeros(shape=(self.nx,self.ny))

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
		
		self.Itot = np.sum(self.stki) / (self.nx * self.ny * obs.nw - num_of_nodes)
		self.Qtot = np.sum(self.stkq) / (self.nx * self.ny * obs.nw - num_of_nodes)
		self.Utot = np.sum(self.stku) / (self.nx * self.ny * obs.nw - num_of_nodes)
		self.Vtot = np.sum(self.stkv) / (self.nx * self.ny * obs.nw - num_of_nodes)

		return



	def read(self, fname, fmt_type=np.float32):
		"""
		Read a binary fortran file

		Parameters
		----------
		fname : str
			Name of the binary file
		fmt_type : type
			Type of the fortran file
			
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
	
	def write(self, fname, fmt_type=np.float32):
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
	
def read_chi2(filename, fmt_type = np.float32):

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
	
