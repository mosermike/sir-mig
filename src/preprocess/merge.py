"""
Merging data downloaded from SDC to a data cube or Hinode data. It puts the 3D arrays into a 4D array to be easier
accessible. In total, four different measurements can be put into 4 data cubes at the same time. For more measurements
in the same folder, the script must be adapted.
"""
import numpy as np
from astropy.io import fits
import os
import sys
sys.path.append("..")
import sir
import definitions as d


def help_page():
	print("merge - Merging data into a data cube")
	print("Usage: python merge.py [OPTION]")
	print()
	sir.option("[1. Pos.]","config file")
	sir.option("[2. Pos.]","Directory of the fits (optional, default ./)")
	sir.option("[3. Pos.]","determines which dataset (001, 002, etc.) (optional, default 001).")
	sys.exit()


def merge(conf, dir, ending):
	"""
	Merges data to a cube

	:param conf: str, config file
	:param dir: str, directory of the fits
	:param ending: str, ending of the dataset (for GRIS data)
	:return: None
	"""
	if conf['instrument'] not in d.instruments:
		print(f"-------> Instrument {conf['instrument']} is not known. Consider adding it to the definitions.py script.")
		print(f"-------> The code assumes that all files in {dir} ending with '.fits' belong to the data.")
		print(f"-------> The code will try to make it work and will ask you for some inputs.")
		print(f"-------> Consider adapting the script yourself or make an issue in the corresponding git website.")

	output_wave = os.path.join(conf['path'], conf['waves'])
	output = os.path.join(conf['path'], conf['cube'])  # where it is saved

	# Files in the folders from different measurements
	filenames = []

	# Add all fits files in a list
	for File in sorted(os.listdir(dir)):
		if File.endswith(".fits"):
			if ("_" + ending + "_0" in File and conf['instrument'] == 'GRIS') or (conf['instrument'] != 'GRIS'):
				filenames.append(dir + File)

	####################
	# Create data cube #
	####################

	example = fits.open(filenames[0])
	stokes = example[0].data
	header = example[0].header

	#########################################
	# 	Load data and store it in an array	#
	#########################################
	# Get no. of wavelengths and pixels ; Reshape the array in a different order
	# ADAPT for another instrument here for the structure
	if conf['instrument'] == 'Hinode':
		nw = stokes.shape[2]
		nx = len(filenames)
		ny = stokes.shape[1]
	elif conf['instrument'] == 'GRIS':
		nw = stokes.shape[1]
		nx = len(filenames)
		ny = stokes.shape[2]
	else:
		temp = input("Specify which entry corresponds to (nw and ny) as a list: [for GRIS: 1,2] ").split("'")
		print("-------> The code assumes that the number of files corresponds to the x direction.")
		nw = stokes.shape[int(temp[0])]
		nx = len(filenames)
		ny = stokes.shape[int(temp[1])]

	print('# Pixels in x  = ', nx)
	print('# Pixels in y  = ', ny)
	print('# Wavelengths  = ', nw)

	data = np.empty((nx, ny, 4, nw))  # shape 4, wavelength, y, 1

	# Create data cube
	for i in range(len(filenames)):
		if i % 100 == 0:
			print(i)

		# Load data from one file
		example = fits.open(filenames[i])

		stokes = example[0].data
		if conf['instrument'] == 'Hinode':
			header = example[0].header
			spbshft = header['SPBSHFT']
			stokes = example[0].data.astype(np.float32)
			stokes = np.transpose(stokes, axes=(1, 0, 2))

			# Correction
			if spbshft > 0:
				stokes[:, 0, :] *= 2
				negative = np.where(stokes[:, 0, :] < 0.0)
				stokes[:, 0, :][negative] += 65536.

			if spbshft > 1:
				stokes[:, 3, :] *= 2
			if spbshft > 2:
				stokes[:, 1, :] *= 2
				stokes[:, 2, :] *= 2

			data[i] = stokes[:, :, :]

		elif conf['instrument'] == 'GRIS':
			stokes = np.transpose(stokes, axes=(3, 2, 0, 1))  # Structure of GRIS stokes, nw, ny, nx
			data[i] = stokes[0, :, :, :]
		else:  # ADAPT for another instrument here for the structure
			print("-------> It assumes the structure as GRIS if it has dim 4.")
			print("-------> It assumes the structure as Hinode if it has dim. 3.")
			if len(stokes.shape) == 4:
				stokes = np.transpose(stokes, axes=(3, 2, 0, 1))  # wanted structure nx, ny, stokes, nw
			elif len(stokes.shape) == 3:
				stokes = np.transpose(stokes, axes=(1, 0, 2))  # wanted structure nx, ny, stokes, nw
			else:
				print("        Manual adaptation necessary!")
				return

			data[i] = stokes[0, :, :, :]

	data = data.reshape(nx, ny, 4, nw)  # shape x, y, values, wavelengths

	# To save data convert to int32
	data = data.astype(np.int32)

	##############################
	# Determine real wavelengths #
	##############################
	if conf['instrument'] == 'GRIS':
		print("[STATUS] Create wavelength arrays assuming the header as for GRIS data")
		ll_a = header["CRVAL3"]
		ll_b = header["CDELT3"]
		llambda = ll_a + ll_b * np.arange(0, header["NAXIS3"])  # Measured wavelength for each pixel
	elif conf['instrument'] == 'Hinode':
		llambda = np.linspace(6300.87730065, 6303.25996109, 112)
	else:
		print(f"Instrument {conf['instrument']} not known. Define wavelength grid manually:")
		mins = input("Lower wavelength in A: ")
		maxs = input("Upper wavelength in A: ")
		values = input("Number of spectral points: ")
		llambda = np.linspace(float(mins), float(maxs), int(values))

	if conf['shift_wave'] != '0':
		print(f"Wavelength Grid shifted by {conf['shift_wave']} mA.")
	llambda = llambda + float(conf['shift_wave']) * 1e-3
	
	np.save(output_wave, llambda)

	print("-------> Saving data (this might take a while) ...")
	if ".npy" in output:
		np.save(output, data)
	else:
		# Write the merged data cube
		hdu = fits.PrimaryHDU(data)
		hdu.header = fits.open(filenames[-1])[0].header
		hdu.writeto(output, overwrite=True)

	print("Saved as \"%s\"" % output)

	# Save some information from the header
	filename = os.path.join(conf['path'], d.header_infos)
	with open(filename, 'w') as f:
		if conf['instrument'] == 'GRIS':
			f.write(f"POINT_ID={header['POINT_ID']}\n")
			f.write(f"DATE-OBS={header['DATE-OBS']}\n")
			f.write(f"CTYPE1={header['CTYPE1']}\n")
			f.write(f"CUNIT1={header['CUNIT1']}\n")
			f.write(f"CRVAL1={header['CRVAL1']}\n")
			f.write(f"CDELT1={header['CDELT1']}\n")
			f.write(f"CUNIT2={header['CUNIT2']}\n")
			f.write(f"CRVAL2={header['CRVAL2']}\n")
			f.write(f"CDELT2={header['CDELT2']}\n")
			f.write(f"CRVAL3={header['CRVAL3']}\n")
			f.write(f"CDELT3={header['CDELT3']}\n")
			f.write(f"NAXIS2={header['NAXIS2']}\n")
			f.write(f"NAXIS3={header['NAXIS3']}\n")
			f.write(f"NAXIS4={header['NAXIS4']}\n")
			f.write(f"POINT_ID={header['POINT_ID']}\n")
			f.write(f"POINT_ID={header['POINT_ID']}\n")
			f.write(f"DSUN_OBS={header['DSUN_OBS']}")
			f.write(f"SHIFT={conf['shift_wave']}\n")
		elif conf['instrument'] == 'Hinode':
			f.write(f"CRVAL1={header['CRVAL1']}\n")
			f.write(f"CDELT1={header['CDELT1']}\n")
			f.write(f"CRVAL2={header['CRVAL2']}\n")
			f.write(f"CDELT2={header['CDELT2']}\n")
			f.write(f"SC_ATTX={header['SC_ATTX']}\n")
			f.write(f"XSCALE={header['XSCALE']}\n")
			f.write(f"XCEN={header['XCEN']}\n")
			f.write(f"YCEN={header['YCEN']}\n")
			f.write(f"CRPIX2={header['CRPIX2']}\n")
			f.write(f"SPBSHFT={header['SPBSHFT']}\n")
			f.write(f"SHIFT={conf['shift_wave']}")
		else:
			print("[NOTE]   No header information file created as instrument is unknown!")

if __name__ == "__main__":

	if "-h" in sys.argv:
		help_page()

	# Read configs and the directory of the data
	conf = sir.read_config(sys.argv[1])
	Dir = sys.argv[2]
	if conf['instrument'] == 'GRIS':
		ends = ''
		if len(sys.argv) > 3 and sys.argv[3][0] != '-':
			ends = sys.argv[3]
		if len(sys.argv) <= 3:
			ends = "001"
		print("-------> The ending of the fits files are", ends)
	else:
		ends = ''

	merge(conf, Dir, ends)
