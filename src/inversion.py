"""
Inversion
=========

Perform inversion of many profiles with simple parallelization using mpi4py.

"""

import os
import shutil
import numpy as np
import time
import sys
from tqdm import tqdm
import glob
from os.path import exists
import datetime

import chi2_stk as c
import create_random_guess as g
import definitions as d
import model_atm as m
import profile_stk as p
import sir




# ***************************************************************************
#									SCATTER									*
#									DATA									*
#****************************************************************************

def scatter_data(conf : dict, comm, rank : int, size : int, return_stk=False):
	"""
	Loads and scatters the data equally into all the processes.

	Parameters
	----------
	conf : Dict
		Config information
	comm : Intercomm
		MPI communicator
	rank : int
		Number of actual process
	size : int
		Number of available processes
	return_stk : bool
		If True, the full stokes profile is also returned

	Returns
	-------
	scatter_data : profile_stk
		Stokes Profiles used in each process
	scatter_data : dict
		Dictionary with the folders, x and y positions for each process
	scatter_data : profile_stk
		Complete stokes profiles in rank 0 (in other ranks 'None') if return_stk = True
	"""
	path = conf["path"]
	
	if rank == 0:
			print("[STATUS] Load and scatter data ...")
			tasks = sir.create_task_folder_list(conf["map"])
			stk1 = p.read_profile(os.path.join(path,conf["cube"]))
			
			if(conf["map"][1] > (stk1.nx-1)) or (conf["map"][3] > (stk1.ny-1)):
				if(conf["map"][1] > (stk1.nx-1)):
					print(f"[scatter_data] Warning: Chosen Map is too big in x! ({conf['map'][1]} > {stk1.nx-1})")
				if(conf["map"][3] > (stk1.ny-1)):
					print(f"[scatter_data] Warning: Chosen Map is too big in y! ({conf['map'][3]} > {stk1.ny-1})")
				print("               Undefined Behaviour! Abort ...")
				sys.exit()


			# Cut data to the wavelength range and to the map
			stk1.cut_to_wave(conf["range_wave"]) # Cut wavelength file to the wished area
			stk1.cut_to_map(conf["map"]) # Cut to the map

			# Create one data cube
			stki = np.copy(stk1.stki).reshape(-1, stk1.nw) # Flatten Stokes I corresponding to the tasks list
			stkq = np.copy(stk1.stkq).reshape(-1, stk1.nw) # Flatten Stokes Q corresponding to the tasks list
			stku = np.copy(stk1.stku).reshape(-1, stk1.nw) # Flatten Stokes U corresponding to the tasks list
			stkv = np.copy(stk1.stkv).reshape(-1, stk1.nw) # Flatten Stokes V corresponding to the tasks list		
			
			waves = np.copy(stk1.wave)

			if not return_stk:
				del stk1 # Free Memory
			
			# Compute the contribution
			stki_size = stki.shape[0]
			chunk_size = stki_size // size
			remainder = stki_size % size # Remainder if number of ranks is not a divisor(?) of the data size
			
			# Divide data into chunks
			stki_chunks = [stki[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
			stkq_chunks = [stkq[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
			stku_chunks = [stku[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
			stkv_chunks = [stkv[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
			tasks_chunks = [tasks['folders'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]
			x_chunks = [tasks['x'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]
			y_chunks = [tasks['y'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]

			# If there is a remainder, distribute the remaining data equally
			if remainder:
				i = 0
				while(i < remainder):
					for j in range(size):
						stki_chunks[j] = np.vstack([stki_chunks[j],stki[-remainder+i]])
						stkq_chunks[j] = np.vstack([stkq_chunks[j],stkq[-remainder+i]])
						stku_chunks[j] = np.vstack([stku_chunks[j],stku[-remainder+i]])
						stkv_chunks[j] = np.vstack([stkv_chunks[j],stkv[-remainder+i]])
									
						tasks_chunks[j].append(tasks['folders'][-remainder+i])
						x_chunks[j] = np.append(x_chunks[j],tasks['x'][-remainder+i])
						y_chunks[j] = np.append(y_chunks[j],tasks['y'][-remainder+i])
						
						i += 1
						if(i >= remainder):
							break
	else:
			stki_chunks = None
			stkq_chunks = None
			stku_chunks = None
			stkv_chunks = None
			x_chunks = None
			y_chunks = None
			tasks_chunks = None
			waves = None
			stk1 = None

	# Scatter data chunks and task structure to all processes
	stki_chunk = comm.scatter(stki_chunks, root=0)
	stkq_chunk = comm.scatter(stkq_chunks, root=0)
	stku_chunk = comm.scatter(stku_chunks, root=0)
	stkv_chunk = comm.scatter(stkv_chunks, root=0)
	task_chunk = comm.scatter(tasks_chunks, root=0)
	x_chunk = comm.scatter(x_chunks, root=0)
	y_chunk = comm.scatter(y_chunks, root=0)
	waves = comm.bcast(waves, root=0)

	# Put the data together
	stk = p.profile_stk(stki_chunk.shape[0],1,nw=stki_chunk.shape[1])
	stk.data_cut_wave = True # Already cut before
	stk._data_cut_map = True # Already cut before
	stk.stki[:,0,:] = stki_chunk
	stk.stkq[:,0,:] = stkq_chunk
	stk.stku[:,0,:] = stku_chunk
	stk.stkv[:,0,:] = stkv_chunk
	stk.wave = waves
	tasks = {
		'folders' : task_chunk,
		'x' : x_chunk,
		'y' : y_chunk
	}
	
	# Free Memory 
	del stki_chunk
	del stkq_chunk
	del stku_chunk
	del stkv_chunk
	del task_chunk
	del x_chunk
	del y_chunk

	if return_stk:
		return stk, tasks, stk1
	else:
		return stk, tasks


def scatter_data_mc(conf : dict, comm, rank : int, size : int, return_stk=False):
	"""
	Scatters the data equally into all the ranks for the MC simulation.

	Parameters
	----------
	conf : Dict
		Config information
	comm : Intercomm
		MPI communicator
	rank : int
		Number of actual process
	size : int
		Number of available processes

	Returns
	-------
	out : Profile
		Stokes Profiles for each process
	out : Dict
		Dictionary with all the folders, x and y position for each process
	"""
	
	if rank == 0:
		print("[STATUS] Load and scatter data ...")
		path = conf["path"]
		if "--no-noise" in sys.argv:
			if rank == 0:
				print("-------> No noise flag used")
				print("-------> Use synthesis profiles")
			stk1 = p.read_profile(os.path.join(path,conf["syn_out"] + d.end_stokes))
		else:
			stk1 = p.read_profile(os.path.join(path,conf["noise_out"] + d.end_stokes))

		tasks = sir.create_task_folder_list(conf["num"])

		if(conf["num"] > (stk1.nx)):
			raise Exception(f"[scatter_data] Warning: Number of Models is too big! ({conf['num']} > {stk1.nx-1})\n               Undefined Behaviour! Abort ...")

		# Create one data cube
		stki = stk1.stki[:,0,:]
		stkq = stk1.stkq[:,0,:]
		stku = stk1.stku[:,0,:]
		stkv = stk1.stkv[:,0,:]
		lines = stk1.indx
		waves = stk1.wave

		if not return_stk:
			del stk1 # Free Memory

		# Compute the contribution
		stki_size = stki.shape[0]
		chunk_size = stki_size // size
		remainder = stki_size % size # Remainder if number of ranks is not a divisor(?) of the data size
		
		# Divide data into chunks
		stki_chunks = [stki[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
		stkq_chunks = [stkq[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
		stku_chunks = [stku[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
		stkv_chunks = [stkv[i*chunk_size:(i+1)*chunk_size] for i in range(size)]
		tasks_chunks = [tasks['folders'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]
		num_chunks = [tasks['x'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]

		# If there is a remainder, distribute the remaining data to the last process
		if remainder:
			i = 0
			while(i < remainder):
				for j in range(size):
					stki_chunks[j] = np.vstack([stki_chunks[j],stki[-remainder+i]])
					stkq_chunks[j] = np.vstack([stkq_chunks[j],stkq[-remainder+i]])
					stku_chunks[j] = np.vstack([stku_chunks[j],stku[-remainder+i]])
					stkv_chunks[j] = np.vstack([stkv_chunks[j],stkv[-remainder+i]])
									
					tasks_chunks[j].append(tasks['folders'][-remainder+i])
					num_chunks[j] = np.append(num_chunks[j],tasks['x'][-remainder+i])
						
					i += 1
					if(i >= remainder):
						break

	else:
		stki_chunks = None
		stkq_chunks = None
		stku_chunks = None
		stkv_chunks = None
		num_chunks = None
		tasks_chunks = None
		waves = None
		lines = None
		stk1 = None

	# Scatter data chunks and task structure to all processes
	stki_chunk = comm.scatter(stki_chunks, root=0)
	stkq_chunk = comm.scatter(stkq_chunks, root=0)
	stku_chunk = comm.scatter(stku_chunks, root=0)
	stkv_chunk = comm.scatter(stkv_chunks, root=0)
	task_chunk = comm.scatter(tasks_chunks, root=0)
	num_chunk = comm.scatter(num_chunks, root=0)
	lines = comm.bcast(lines, root=0)
	waves = comm.bcast(waves, root=0)

	# Put the data together
	stk = p.profile_stk(stki_chunk.shape[0],1,stki_chunk.shape[1])
	stk.data_cut = True # Already cut before
	stk.indx = lines
	stk.stki[:,0] = stki_chunk
	stk.stkq[:,0] = stkq_chunk
	stk.stku[:,0] = stku_chunk
	stk.stkv[:,0] = stkv_chunk
	stk.wave = waves
	tasks = {
		'folders' : task_chunk,
		'num' : num_chunk,
	}
	
	# Free Memory 
	del stki_chunk
	del stkq_chunk
	del stku_chunk
	del stkv_chunk
	del task_chunk
	del num_chunk

	if return_stk:
		return stk, tasks, stk1
	else:
		return stk, tasks

"""
*****************************************************************************
	*							EXECUTION OF								*
*								INVERSIONS									*
*****************************************************************************
"""

def execute_inversion_1c(conf : dict, task_folder : dict, rank : int):
	"""
	Executes inversion and creates if needed random guesses. Make sure that
	when this function is entered, the os is in the right directory! This
	function can be used for the mode 'MC' and '1C'.
	
	Parameters
	----------
	config : dict
		Dictionary containing the configuration of the simulation
	task_folder : str
		String of the directory where this code should executed
	rank : int
		Process Number
	
	Returns
	-------
	out : float
		Best chi2 value of the inversion(s)
	
	"""
	# Change the directory
	old_dir = os.path.abspath(os.getcwd())
	os.chdir(task_folder)

	# Define parameters for simplicity reasons
	shutil.copy(conf['model'], d.model_inv)
	model = d.model_inv
	guess1 = d.model_inv.replace(".mod", "")  # For simplicity
	cycles = conf["cycles"]
	chi_file = d.inv_trol_file[:d.inv_trol_file.rfind('.')] + ".chi"

	#########################
	# START INVERSION #
	#########################
	# Determine whether random guesses are wished or not
	if conf["random_guess"] > 0:
		# Initialize the chi2 map
		chi2 = np.zeros(conf["random_guess"])

		# Perform inversion for each guess model and copy files
		for i in range(conf["random_guess"]):
			it = 0
			while it < 50: # stop (repeating) inversion when chi2 file is correctly generated => inversion finished correctly or more than 50 iterations
				# Create New Guess
				g.create_guesses_1c(conf, output="./", number=i+1)
				# Copy to the model
				shutil.copy(d.model + str(i+1) + ".mod", d.model_inv)
				# Execute inversion
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Sometimes it does not converge and when only 8 cols are used, the programme could crash
				# => results into NaNs
				# Also other bugs/problems (SVD singularity) could result into failure => Check dimension of inv.chi file
				if exists(chi_file):  # Chi file exists
					if np.genfromtxt('inv.chi').shape == (int(cycles), 2) or int(cycles) == 1:  # All the cycles finished
						break  # Get out of the while loop
					else:  # Do inversion again
						os.remove(chi_file)
				else: # No inv.chi file generated => There might be a problem with sir.x => Break loop after it is greater than 50
					it += 1
			# If it is greater than 50, there might be something wrong with sir.x
			if it >= 50:
				# Print last log entry
				if rank == 0:
					print(f"[ERROR] Check your sir.x file and the log file in the {d.task_start} folders. There might be a problem with sir.x")
					with open('inv.log') as f:
						lines = []
						for line in f:
							lines.append(line)
						print("[LAST LOG ENTRIES]: ")
						print("------------------------------ ")
						print(lines[-5])
						print(lines[-4])
						print(lines[-3])
						print(lines[-2])
						print(lines[-1])
						print("------------------------------ ")
					del lines

					sys.exit()
			# Read the chi2 file
			chi = sir.read_chi2(chi_file, task_folder)

			# Inversion did not work => Put to a significantly high value to not affect the choice
			# of the best inversion if another converged
			if chi == 0.0:
				chi = 1e9

			chi2[i] = chi # Store the chi2 value

			# Move files to make space for new inversions
			shutil.move(chi_file, chi_file + "__" + str(i+1))
			shutil.move(f"{guess1}_{cycles}.mod", f"{guess1}__{str(i+1)}.mod")
			shutil.move(f"{guess1}_{cycles}.err", f"{guess1}__{str(i+1)}.err")
			shutil.move(f"{guess1}_{cycles}.per", f"{guess1}__{str(i+1)}.per")

			# Remove not needed files to make space for the next inversion
			for j in range(int(cycles)-1):
				os.remove(f"{guess1}_{j+1}.mod")
				os.remove(f"{guess1}_{j+1}.err")
				os.remove(f"{guess1}_{j+1}.per")

		# Check if the best model has NaN entries
		index = chi2.argmin()+1
		mod_is_nan = np.isnan(np.loadtxt(f"{guess1}__{str(index)}.mod", skiprows=1).transpose()).any()
		
		# All inversions failed => Redo inversion until one does not fail
		if chi2.min() > 1e8 or mod_is_nan:
			it = 0 # Number of iterations until sir converges
			chi2_best = 0.0
			while chi2_best == 0.0 or mod_is_nan:
				# Copy the guess which did not work
				shutil.copy(model, "bad_guess" + str(i+1) + ".mod")

				# Create New Guess
				g.create_guesses_1c(conf, output="./", number=i+1)

				# Copy to the model
				shutil.copy(d.model + str(i+1) + ".mod", d.model_inv)

				# Execute inversion again
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Read chi2 value and the model
				chi2_best = sir.read_chi2(chi_file, task_folder)
				mod_is_nan = np.isnan(np.loadtxt(f"{guess1}_{cycles}.mod", skiprows=1)).any()

				# +1 iteration
				it += 1

				# break the iteration
				if it > 50:
					print(task_folder[task_folder.rfind('/')+1:], " failed!")
					break
			
			shutil.copy(f"{d.model_inv}", f"{d.best_guess}")  # Copy best guess model
			shutil.move(f"{guess1}_{cycles}.mod", f"best.mod")
			shutil.move(f"{guess1}_{cycles}.err", f"best.err")
			shutil.move(f"{guess1}_{cycles}.per", f"best.per")
			# Warn if the repetition of the inverse needed to be done more than 20 times
			if it > 10:
				print(f"\nNote that the while loop run in {task_folder[task_folder.rfind('/')+1:]} took {it} iterations because chi2 was 0.0 or the model contained NaNs...")

		# Determine best fit and copy stuff
		else:

			chi2_best = chi2.argmin()
			shutil.copy(f"inv.chi__{str(chi2_best+1)}",chi_file)
			shutil.move(f"{guess1}__{str(chi2_best+1)}.mod", f"best.mod")
			shutil.move(f"{guess1}__{str(chi2_best+1)}.err", f"best.err")
			shutil.move(f"{guess1}__{str(chi2_best+1)}.per", f"best.per")
			shutil.copy(f"{d.model}{str(chi2_best+1)}.mod", f"{d.best_guess}") # Copy best guess model

	else:
		# Check if vmacro is different in the header of the model and if yes change it
		header = np.loadtxt(conf['model'], max_rows=1)
		if header[0] != round(float(conf['vmacro']),4):
			temp_mod = m.read_model(conf['model'])
			temp_mod.vmacro[0,0] = float(conf['vmacro'])
			temp_mod.write_model(conf['model'], 0, 0)
		
		# Copy to inversion model
		shutil.copy(conf["model"], d.model_inv)

		# Perform inversion once and read chi2 value
		os.system(f"echo {d.inv_trol_file} | ./sir.x > inv.log")
		
		# Check if chi2 file exists and inversion finished successfully
		if not exists(chi_file):
			print("[read_chi2] " + chi_file + " does not exist in " + task_folder + ".")
			print(f"[ERROR] Check your sir.x file and the log file in the {d.task_start} folders. There might be a problem with sir.x")
			with open('inv.log') as f:
				lines = []
				for line in f:
					lines.append(line)
				print("[LAST LOG ENTRIES]: ")
				print("------------------------------ ")
				print(lines[-5])
				print(lines[-4])
				print(lines[-3])
				print(lines[-2])
				print(lines[-1])
				print("------------------------------ ")
			del lines
			sys.exit(1)
			

		shutil.move(f"{d.guess.replace('.mod','')}_{cycles}.mod", f"best.mod")
		shutil.move(f"{d.guess.replace('.mod','')}_{cycles}.err", f"best.err")
		shutil.move(f"{d.guess.replace('.mod','')}_{cycles}.per", f"best.per")
		shutil.move(d.model_inv, d.best_guess)
		# Remove not needed files to make space for the next inversion
		for j in range(int(cycles)-1):
			os.remove(f"{guess1}_{j+1}.mod")
			os.remove(f"{guess1}_{j+1}.err")
			os.remove(f"{guess1}_{j+1}.per")

	# Change the directory
	os.chdir(old_dir)
	return


def execute_inversion_2c(conf : dict, task_folder : dict, rank : int):
	"""
	Executes inversion and creates if needed random guesses. Make sure that
	when this function is entered, the os is in the right directory! This is
	used for the mode '2C'.
	
	Parameters
	----------
	conf : Dictionary
		Dictionary containing the configuration of the simulation
	task_folder : str
		Task folder where the code is going
	rank : int
		Process Number

	Returns
	-------
	None

	Raise
	-----
	FileNotFoundError
		if a file is not found
	"""
	old_pw = os.path.abspath(os.getcwd())
	os.chdir(task_folder)

	# Define parameters for simplicity reasons
	model1 = conf['model1'] + "1" # Adding 1 in case models have the same name
	model2 = conf['model2'] + "2" # Adding 2 in case models have the same name
	cycles = conf["cycles"]
	chi_file = d.inv_trol_file[:d.inv_trol_file.rfind('.')] + ".chi"

	####################################
	#		START INVERSION		#
	####################################
	# Determine whether random guesses are wished or not
	if conf["random_guess"] > 0:
		# Initialize the chi2 map
		chi2 = np.zeros(conf["random_guess"])

		# Perform inversion for each guess model and copy files
		for i in range(conf["random_guess"]):
			it = 0
			while it < 50: # Only stops an inversion, if the chi file is correct => inversion finished properly
				# Create New Guess
				g.create_guesses_2c(conf, output = "./", number = i+1)	
				if not exists(f'{d.model1}' + str(i+1) + ".mod"):
					raise FileNotFoundError(f'No model {i+1} in {task_folder} with the name {d.model1}{i+1}.mod')
				if not exists(f'{d.model2}' + str(i+1) + ".mod"):
					raise FileNotFoundError(f'No model {i+1} in {task_folder} with the name {d.model2}{i+1}.mod')
				# Copy to the model
				shutil.copy(d.model1 + str(i+1) + ".mod",d.guess1)
				shutil.copy(d.model2 + str(i+1) + ".mod",d.guess2)
				# Execute inversion
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Sometimes it does not converge and when only 8 cols are used, the programm could crash => results into NaNs
				# Also other bugs/problems (SVD singularity) could result into failure => Check dimension of inv.chi file
				if exists(chi_file): # Chi file exists
					if  np.genfromtxt('inv.chi').shape == (int(cycles),2) or int(cycles) == 1: # All the cycles finished
						break # Get out of the while loop
					else: # Do inversion again
						os.remove(chi_file)
				else: # No inv.chi file generated => There might be a problem with sir.x => Break loop after it is greater than 50
					it += 1
			# If it is greater than 50, there might be something wrong with sir.x
			if it >= 50:
				if rank == 0:
					print(f"[ERROR] Check your sir.x file and the log file in the {d.task_start} folders. There might be a problem with sir.x")
					# Print last log entry
					with open('inv.log') as f:
						lines = []
						for line in f:
							lines.append(line)
						print("[LAST LOG ENTRIES]: ")
						print("------------------------------ ")
						print(lines[-5])
						print(lines[-4])
						print(lines[-3])
						print(lines[-2])
						print(lines[-1])
						print("------------------------------ ")
						del lines
					sys.exit()

			# Read chi2 file
			chi = sir.read_chi2(chi_file, task=task_folder)

			# Inversion did not work if chi2 is 0.0 => Put to a significantly high value to not affect the choice
			# of the best inversion if another converged
			if chi == 0.0:
				chi = 1e9

			chi2[i] = chi # Store the chi2 value

			# Move files to make space for new inversions
			shutil.move(chi_file, chi_file + "__" + str(i+1))
			shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.mod", f"{d.model1}_{str(i+1)}.mod")
			shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.err", f"{d.model1}_{str(i+1)}.err")
			shutil.move(f"{d.guess2.replace('.mod','')}_{cycles}.mod", f"{d.model2}_{str(i+1)}.mod")
			shutil.move(f"{d.guess2.replace('.mod','')}_{cycles}.err", f"{d.model2}_{str(i+1)}.err")
			shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.per", f"{d.model1}_{str(i+1)}.per")

			# This for loop might be removed in a future release to save computing time
			for i in range(int(cycles)-1):
				os.remove(f"{d.guess1.replace('.mod','')}_{i+1}.mod")
				os.remove(f"{d.guess1.replace('.mod','')}_{i+1}.err")
				os.remove(f"{d.guess1.replace('.mod','')}_{i+1}.per")
				os.remove(f"{d.guess2.replace('.mod','')}_{i+1}.mod")
				os.remove(f"{d.guess2.replace('.mod','')}_{i+1}.err")

		# Check if the model has NaNs
		index = str(chi2.argmin()+1)
		mod_is_nan1 = np.isnan(np.loadtxt(f"{d.model1}_{index}.mod", skiprows=1)).any()
		mod_is_nan2 = np.isnan(np.loadtxt(f"{d.model2}_{index}.mod", skiprows=1)).any()

		# All inversions failed => Redo inversion until one does not fail
		if chi2.min() > 1e8 or mod_is_nan1 or mod_is_nan2:
			it = 0 # Number of iterations until sir converges
			chi2_best = 0.0
			while chi2_best == 0.0 or mod_is_nan1 or mod_is_nan2:
				# Create New Guess
				g.create_guesses_2c(conf, output = "./", number = 1)
				# Copy to the model
				shutil.copy(f"{d.model1}1.mod",d.guess1)
				shutil.copy(f"{d.model2}1.mod",d.guess2)

				# Execute inversion again
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Read chi2 value and model
				chi2_best = sir.read_chi2(chi_file, task=task_folder)
				mod_is_nan1 = np.isnan(np.loadtxt(f"{d.guess1.replace('.mod','')}_{cycles}.mod", skiprows=1)).any()
				mod_is_nan2 = np.isnan(np.loadtxt(f"{d.guess2.replace('.mod','')}_{cycles}.mod", skiprows=1)).any()

				# +1 iteration
				it += 1

				if (it > 50):
					print(task_folder[task_folder.rfind('/')+1:], " failed!")
					break
			
			shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.mod", f"best1.mod") # Model for Model 1
			shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.err", f"best1.err") # Error for Model 1
			shutil.move(f"{d.guess2.replace('.mod','')}_{cycles}.mod", f"best2.mod") # Model for Model 2
			shutil.move(f"{d.guess2.replace('.mod','')}_{cycles}.err", f"best2.err") # Error for Model 2
			shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.per", f"best.per")  # Profile
			shutil.move(f"{d.guess1}", f"{d.best_guess1}") # Copy best guess model 1
			shutil.move(f"{d.guess2}", f"{d.best_guess2}") # Copy best guess model 2

			# Warn if the repetition of the inverse needed to be done more than 20 times
			if it > 10:
				print(f"\nNote that the while loop run in {task_folder[task_folder.rfind('/')+1:]} took {it} iterations because chi2 was 0.0...")

		# Determine best fit and copy stuff
		else:

			chi2_best = chi2.argmin()
			shutil.move(f"inv.chi__{str(chi2_best+1)}",chi_file)
			shutil.move(f"{d.model1}_{chi2_best+1}.mod", f"best1.mod") # Model for Model 1
			shutil.move(f"{d.model1}_{chi2_best+1}.err", f"best1.err") # Error for Model 1
			shutil.move(f"{d.model2}_{chi2_best+1}.mod", f"best2.mod") # Model for Model 2
			shutil.move(f"{d.model2}_{chi2_best+1}.err", f"best2.err") # Error for Model 2
			shutil.move(f"{d.model1}_{chi2_best+1}.per", f"best.per")  # Profile
			shutil.move(f"{d.model1}{chi2_best+1}.mod", f"{d.best_guess1}") # Copy best guess model 1
			shutil.move(f"{d.model2}{chi2_best+1}.mod", f"{d.best_guess2}") # Copy best guess model 2
			chi2_best = chi2[chi2_best]

	else:
		# Correct for vmacro
		header = np.loadtxt(model1, max_rows=1)
		if header[0] != round(float(conf['vmacro']),4) or header[1] != round(float(conf['fill'].split(",")[0]),4):
			temp_mod = m.read_model(model1)
			temp_mod.vmacro[0,0] = float(conf['vmacro'])
			temp_mod.fill[0,0] = float(conf['fill'].split(",")[0])
			temp_mod.write_model(model1, 0, 0)
		
		# Correct of different vmacro
		header = np.loadtxt(model2, max_rows=1)
		if header[0] != round(float(conf['vmacro']),4) or header[1] != round(float(conf['fill'].split(",")[1]),4):
			temp_mod = m.read_model(model2)
			temp_mod.vmacro[0,0] = float(conf['vmacro'])
			temp_mod.fill[0,0] = float(conf['fill'].split(",")[1])
			temp_mod.write_model(model2, 0, 0)

		# Perform inversion once and read chi2 value
		shutil.copy(model1,d.guess1)
		shutil.copy(model2,d.guess2)
		os.system(f"echo {d.inv_trol_file} | ./sir.x > inv.log")
		# Check if inversion was successful
		if not exists(chi_file):
			print("[read_chi2] " + chi_file + " does not exist in " + task_folder + ".")
			print(f"[ERROR] Check your sir.x file and the log file in the {d.task_start} folders. There might be a problem with sir.x")
			with open('inv.log') as f:
				lines = []
				for line in f:
					lines.append(line)
				print("[LAST LOG ENTRIES]: ")
				print("------------------------------ ")
				print(lines[-5])
				print(lines[-4])
				print(lines[-3])
				print(lines[-2])
				print(lines[-1])
				print("------------------------------ ")
			del lines
			sys.exit(1)
		
		shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.mod", f"best1.mod")
		shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.err", f"best1.err")
		shutil.move(f"{d.guess2.replace('.mod','')}_{cycles}.mod", f"best2.mod")
		shutil.move(f"{d.guess2.replace('.mod','')}_{cycles}.err", f"best2.err")
		shutil.move(f"{d.guess1.replace('.mod','')}_{cycles}.per", f"best.per")
		shutil.copy(f"{model1}", f"{d.best_guess1}") # Copy best guess model 1
		shutil.copy(f"{model2}", f"{d.best_guess2}") # Copy best guess model 2
		
		# This for loop might be removed in a future release to save computing time => model.py needs to be adjusted then
		for i in range(int(cycles)-1):
			os.remove(f"{d.guess1.replace('.mod','')}_{i+1}.mod")
			os.remove(f"{d.guess1.replace('.mod','')}_{i+1}.err")
			os.remove(f"{d.guess1.replace('.mod','')}_{i+1}.per")
			os.remove(f"{d.guess2.replace('.mod','')}_{i+1}.mod")
			os.remove(f"{d.guess2.replace('.mod','')}_{i+1}.err")

	# Change to old path
	os.chdir(old_pw)

	return 



"""
*****************************************************************************
*							ONE COMPONENT INVERSION							*
*									START									*
*****************************************************************************
"""


def inversion_1c(conf : dict, comm, rank : int, size : int, MPI, debug=False, progress=True):
	"""
	Performs the inversion of all the models for the 1 component inversion.

	Parameters
	----------
	config : dict
		Dictionary with all information from the config file
	comm : Intercomm
		MPI communicator
	rank : int
		Number of actual process
	size : int
		Number of available processes
	MPI : library
		 MPI library imported as 'from mpi4py import MPI'
	debug : bool,optional
		Execute in debugging mode
	progress : bool,optional
		Print a progress bar

	Returns
	-------
	None
	"""
	start = time.time()
	####################################
	# READ PARAMETERS FROM CONFIG FILE #
	####################################
	# Define parameters for easier understanding
	path = conf["path"]
	model = conf["model"]

	line_file = conf['line'] # line file

	# Read inversion stuff
	Map = conf['map']
	if rank == 0:
		print(f"-------> Selected Map ({Map[0]},{Map[1]},{Map[2]},{Map[3]})")

	abundance_file = conf['abundance']  # Abundance file
	
	# Write psf function, if needed
	if rank == 0:
		if conf['psf'] != '':
			#print("-------> Spectral PSF is used")
			if conf['psf'].split(" ")[0] == "gauss":
				print(r"-------> Create Gaussian spectral PSF with σ = " + conf['psf'].replace("gauss ","") + " mÅ")
				sir.write_gauss_psf(float(conf['psf'].replace("gauss ","")), os.path.join(path,d.psf))
			else:
				print(r"-------> Use spectral PSF file {conf['psf']}")
				if(conf['psf'] != d.psf):
					shutil.copy(os.path.join(conf['path'],conf['psf']),os.path.join(conf['path'],d.psf))
	
	# Create guess from bin file if wanted
	if conf["guess"] != '':
		guess = None
		if rank == 0:
			print(f"-------> File {conf['guess'].split(' ')[0]} used as initial guess/base model")
			guess = m.read_model(os.path.join(path, conf["guess"].split(" ")[0]))
			if guess.nx > Map[1]-Map[0]+1 or guess.ny > Map[3]-Map[2]+1:
					print("[Warn]   The shapes of the initial guess/base model are bigger than the map in the config.")
			if guess.nx < Map[1]-Map[0]+1 or guess.ny < Map[3]-Map[2]+1:
					print(f"[ERROR]  The shapes ({guess.nx},{guess.ny}) of the initial guess/base model are not compatible with the map ({Map[1]-Map[0]+1},{Map[3]-Map[2]+1}) in the config.")
					return
			
			# Small variation around the guess
			if(len(conf["guess"].split(" ")) > 1):
				factor = float(conf["guess"].split(" ")[1])
				print(f"-------> Create small variation around the initial guess with the factor {factor}")
				guess = g.create_small_guess(guess, conf["random_pars"], factor)
		
		guess = comm.bcast(guess, root=0)
		
	if rank == 0:
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)
	
	if rank == 0:
		# Write which parameters are randomised
		if conf["random_guess"] > 0 or (conf["guess"] != '' and len(conf["guess"].split(" ")) > 1):
			random_pars = conf["random_pars"]
			print("-------> Parameters to be randomised: ", end='')
			print_out = ''
			if 'T' in random_pars:
				print_out += "Temperature, "
			if 'Pe' in random_pars:
				print_out += "El. Pressure, "
			if 'vmicro' in random_pars:
				print_out += "Microturbulent Velocity, "
			if 'B' in random_pars:
				print_out += "Magnetic Field, "
			if 'vlos' in random_pars:
				print_out += "Line-of-Sight Velocity, "
			if 'gamma' in random_pars:
				print_out += "Inclination, "
			if 'phi' in random_pars:
				print_out += "Azimuth, "
			if 'z' in random_pars:
				print_out += "Height, "
			if 'Pg' in random_pars:
				print_out += "Gas Pressure, "
			if 'rho' in random_pars:
				print_out += "Density, "
			print(print_out[0:-2] + ".")

			if conf['random_guess'] == 1:
				print(f"-------> Create {conf['random_guess']} random guess per pixel.")
			else:
				print(f"-------> Create {conf['random_guess']} random guesses per pixel.")
		elif conf["guess"] == '':
			print(f"-------> Use Base Model '{conf['model']}' as initial guess.")

	########################
	# START INVERSION PART #
	########################
	total_jobs = 1  # Total performed jobs across all processes
	tasks = sir.create_task_folder_list(Map) # Structure tasks
	max_jobs = len(tasks['folders'])  # For comm.allreduce function

	#########################
	# SCATTER AND LOAD DATA #
	#########################
	# Load and scatter data => Saving memory and time
	if conf['chi2'] == "1":
		stk, tasks, obs = scatter_data(conf, comm, rank, size, True)
	else:
		stk, tasks = scatter_data(conf, comm, rank, size, False)

	comm.barrier()
	
	# Write the control file with the information from the config file
	if rank == 0:
		print("-------> Write control and grid file")
		sir.write_control(os.path.join(conf['path'],d.inv_trol_file), conf)
		# Write Grid file based on the chosen wavelength ranges in the config file
		sir.write_grid(conf, os.path.join(path,d.Grid), stk.wave)

	if rank == 0:
		print("[STATUS] Start Computing Inversions ...")

	# Track local progress
	performed_models = 0  # Counts how many models are performed

	# Root process initializes the progress bar
	if rank == 0 and progress:
		pbar = tqdm(total=max_jobs, desc="Overall Progress", file=sys.stdout)
	comm.barrier()

	for i in range(0, len(tasks['folders'])):
		####################################
		# Create the folder and copy stuff #
		####################################
		# Read x and y position, create task folder
		task_folder = os.path.join(path, tasks['folders'][i])
		x = tasks['x'][i]
		y = tasks['y'][i]
		os.makedirs(task_folder, exist_ok=True)
					
		# Write the data from the cube into a profile file for SIR
		stk.write_profile(os.path.join(path, task_folder, d.profile_obs), i, 0, os.path.join(conf['path'],d.Grid))
				
		# Copy stuff for the inversion
		if conf['psf'] != '':
			sir_files = [d.inv_trol_file, model, "sir.x", line_file, d.Grid, abundance_file, d.psf]
		else:
			sir_files = [d.inv_trol_file, model, "sir.x", line_file, d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(task_folder, sir_file))

		# Create guess from npy file if wanted
		# Copy the data from the imported numpy array to the model in the task folder
		if conf["guess"] != '':
			# Write the new initial model from here:
			x1 = x - Map[0]  # x position taking the Map into account
			y1 = y - Map[2]  # y position taking the Map into account
			guess.write_model(os.path.join(task_folder, model), x1, y1)

		#####################
		# Perform inversion #
		#####################
		# Create random guesses and select best value and perform inversion
		execute_inversion_1c(conf, task_folder, rank)
		
		# Update performed_models
		performed_models += 1

		# Do not do allreduce for the last step as the code does not move on from here
		if total_jobs < (max_jobs - max_jobs % size):
			# Root process updates the progress bar
			total_jobs = comm.allreduce(performed_models, op=MPI.SUM)
		
		# Update progres bar
		if rank == 0 and progress:
			pbar.n = total_jobs
			pbar.refresh()
		elif rank == 0:
			print(f"\rTotal Jobs: {total_jobs}/{max_jobs}")
	
	comm.barrier()
		
	if rank == 0 and progress:
		pbar.n = max_jobs
		pbar.refresh()
		pbar.close()
	elif rank == 0:
			print(f"\rTotal Jobs: {total_jobs}/{max_jobs}")
			
	wave = np.copy(stk.wave)
	del stk # Free Memory

	comm.barrier()
		
	##################################################
	# Read all the results and put it into npy files #
	##################################################
	if rank == 0:
		print("[STATUS] Gathering results...")

		# Redefine tasks as now all the tasks are read
		tasks = sir.create_task_folder_list(Map) # Structure tasks

		# Create directory if / exists in inv_out in config
		if '/' in conf['inv_out']:
			temp = os.path.join(path, conf["inv_out"])
			if not exists(temp[:temp.rfind('/')]):
				os.mkdir(temp[:temp.rfind('/')])		
		

		print("-------> Read and Write Profiles ...")
		stokes_inv = p.profile_stk(0,0,0)
		stokes_inv.wave = wave # Copy wavelength positions
		stokes_inv.read_results(tasks, "best.per", path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		stokes_inv.write(os.path.join(path,conf['inv_out']) + d.end_stokes)

		if conf["chi2"] == "":
			del stokes_inv
			
		print("-------> Read and Write Models ...")
		models_inv = m.model_atm(0,0,0)
		errors_inv = m.model_atm(0,0,0)
		best_guesses = m.model_atm(0,0,0)

		models_inv.read_results(tasks, 'best.mod', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		errors_inv.read_results(tasks, 'best.err', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		best_guesses.read_results(tasks, d.best_guess, path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)

		models_inv.write(os.path.join(path,conf['inv_out']) + d.end_models)
		errors_inv.write(os.path.join(path,conf['inv_out']) + d.end_errors)
		best_guesses.write(os.path.join(path,conf['inv_out'] + d.best_guess_file))

		del models_inv
		del errors_inv
		del best_guesses

		if conf['chi2'] == "1":
			chi2 = c.chi2_stk(stokes_inv.nx,stokes_inv.ny)
			print("-------> Compute χ² (this might take a while) ...")
			
			# Number of Nodes in the last step
			num_of_nodes = 0
			if conf['nodes_temp'] != "":
				num_of_nodes+= int(conf['nodes_temp'].split(",")[-1])
			if conf['nodes_magn'] != "":
				num_of_nodes+= int(conf['nodes_magn'].split(",")[-1])
			if conf['nodes_vlos'] != "":
				num_of_nodes+= int(conf['nodes_vlos'].split(",")[-1])
			if conf['nodes_gamma'] != "":
				num_of_nodes+= int(conf['nodes_gamma'].split(",")[-1])
			if conf['nodes_phi'] != "":
				num_of_nodes+= int(conf['nodes_phi'].split(",")[-1])
			
			# Compute chi2
			chi2.compute(obs, stokes_inv, [float(i) for i in conf["weights"]], num_of_nodes)
			print("-------> Total χ² = %.3f" % chi2.total)
			print("-------> Write χ² ...")
			chi2.write(os.path.join(path, conf['inv_out'] + d.end_chi2))
			del obs
			del stokes_inv
			del num_of_nodes
			del chi2

		if not debug:
			print("-------> Delete created SIR files")
			# Delete the folder
			for i in range(len(tasks['x'])):
				# Remove folder
				shutil.rmtree(os.path.join(path,tasks['folders'][i]))

			if conf['psf'] != '' and conf['psf'] != d.psf:
				os.remove(os.path.join(path,d.psf))

			os.remove(os.path.join(path,d.inv_trol_file))
			os.remove(os.path.join(path,d.Grid))
		else:
			print("-------> No created files are deleted. Debug option active.")
		
		# Print needed time
		end = time.time()
		Time = str(datetime.timedelta(seconds=end-start)).split(".")[0]
		print(f"-------> Finished in {Time}")

	comm.barrier()

	return



"""
*****************************************************************************
*							MONTE CARLO SIMULATION							*
*									START									*
*****************************************************************************
"""

def inversion_mc(conf : dict, comm, rank : int, size : int, MPI, debug=False,progress=True):
	"""
	Performs the inversion of all the models for the MC simulation.

	Parameters
	----------
	conf : dict
		Dictionary with all information from the config file
	comm : Intercomm
		MPI communicator
	rank : int
		Number of actual process
	size : int
		Number of available processes
	MPI : library
		MPI library imported as 'from mpi4py import MPI'
	debug : bool,optional
		Execute in debugging mode
	progress : bool,optional
		Print a progress bar

	Returns
	-------
	None

	"""
	start = time.time()
	####################################
	# READ PARAMETERS FROM CONFIG FILE #
	####################################
	# Define parameters for easier understanding
	path = conf["path"]
	model = conf["model"]
	line_file = conf['line']  # line file
	
	# Write the control file with the information from the config file
	if rank == 0:
		print(f"-------> Write control file")
		sir._write_control_mc(os.path.join(conf['path'], d.inv_trol_file), conf, "inv")
	abundance_file = conf['abundance']  # Abundance file
	

	# Create guess from bin file if wanted
	if conf["guess"] != '':
		guess = None
		if rank == 0:
			guess = m.read_model(os.path.join(path, conf["guess"].split(" ")[0]))
			print(f"-------> File {conf['guess'].split(' ')[0]} used as initial guess/base model")
			if guess.nx > conf['num']:
					print("[Warn] The shapes of the initial guess/base model are bigger than the map in the config.")
			if guess.nx < conf['num']:
					print(f"[ERROR]  The shapes ({guess.nx},{guess.ny}) of the initial guess/base model are not compatible with the map ({conf['num']},1) in the config.")
					return
			
			# Small variation around the guess
			if(len(conf["guess"].split(" ")) > 1):
				factor1 = float(conf["guess"].split(" ")[1])
				print(f"-------> Create small variation around the initial guess with the factor {factor1}")
				guess = g.create_small_guess(guess, conf["random_pars"], factor1)
		
		guess = comm.bcast(guess, root=0)

	####################
	# CREATE GRID FILE #
	####################
	if rank == 0:
		# Write Grid file based on the chosen wavelength ranges in the config file
		sir.write_grid(conf, os.path.join(path, d.Grid))
		# Write which parameters are randomised
		if conf["random_guess"] > 0 or (conf["guess"] != '' and len(conf["guess"].split(" ")) > 1):
			random_pars = conf["random_pars"]
			print("-------> Parameters to be randomised: ", end='')
			print_out = ''
			if 'T' in random_pars:
				print_out += "Temperature, "
			if 'Pe' in random_pars:
				print_out += "El. Pressure, "
			if 'vmicro' in random_pars:
				print_out += "Microturbulent Velocity, "
			if 'B' in random_pars:
				print_out += "Magnetic Field, "
			if 'vlos' in random_pars:
				print_out += "Line-of-Sight Velocity, "
			if 'gamma' in random_pars:
				print_out += "Inclination, "
			if 'phi' in random_pars:
				print_out += "Azimuth, "
			if 'z' in random_pars:
				print_out += "Height, "
			if 'Pg' in random_pars:
				print_out += "Gas Pressure, "
			if 'rho' in random_pars:
				print_out += "Density, "
			print(print_out[0:-2] + ".")
			if conf['random_guess'] == 1:
				print(f"-------> Create {conf['random_guess']} random guess per pixel.")
			else:
				print(f"-------> Create {conf['random_guess']} random guesses per pixel.")
		else:
			if rank == 0:
				print(f"-------> Use Base Model '{conf['model']}' as initial guess.")

	if rank == 0:
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)

	comm.barrier()

	########################
	# START INVERSION PART #
	########################
	if rank == 0:
		print("[STATUS] Start Inversions")
	performed_models = 0  # Counts how many models are performed
	total_jobs = 1  # Total performed jobs across all processes
	max_jobs = conf['num']  # For comm.allreduce function

	if conf['chi2'] == "1":
		stk, tasks, obs = scatter_data_mc(conf, comm, rank, size, True)
	else:
		stk, tasks = scatter_data_mc(conf, comm, rank, size, False)

	# Root process initializes the progress bar
	if rank == 0 and progress:
		pbar = tqdm(total=max_jobs, desc="Overall Progress", file=sys.stdout)
		
	comm.barrier()

	for i in range(0,len(tasks['num'])):
		####################################
		# Create the folder and copy stuff #
		####################################
		# Create task folder
		task_folder = os.path.join(path, tasks['folders'][i])
		os.makedirs(task_folder, exist_ok=True)

		# Write the data from the cube into a profile file for SIR
		stk.write_profile_mc(os.path.join(task_folder, d.profile), i)

		sir_files = [d.inv_trol_file, model, "sir.x", line_file, d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(task_folder, sir_file))
		
		# Create guess from bin file
		if conf["guess"] != '':
			# Write the new initial model from here:
			guess.write_model(os.path.join(task_folder, model), i, 0)

		###############################
		# 	Perform inversion		#
		###############################
		# Create random guesses and select best value and perform inversion
		execute_inversion_1c(conf, task_folder, rank)
		
		performed_models += 1
	
		# Do not do allreduce for the last step as the code does not move on from here
		if total_jobs < (max_jobs - max_jobs % size):
			# Root process updates the progress bar
			total_jobs = comm.allreduce(performed_models, op=MPI.SUM)
		
		# Update progres bar
		if rank == 0 and progress:
			pbar.n = total_jobs
			pbar.refresh()
		elif rank == 0:
			print(f"\rTotal Jobs: {total_jobs}/{max_jobs}", end="")

	comm.barrier()
		
	if rank == 0 and progress:
		pbar.n = max_jobs
		pbar.refresh()
		pbar.close()
	elif rank == 0:
			print(f"\rTotal Jobs: {total_jobs}/{max_jobs}")
	##################################################
	# Read all the results and put it into npy files #
	##################################################
	if rank == 0:
		print("[STATUS] Gathering results...")

		tasks = sir.create_task_folder_list(conf["num"])

		# Create directory if / exists in inv_out in config
		if '/' in conf['inv_out']:
			temp = os.path.join(path, conf["inv_out"])
			if not exists(temp[:temp.rfind('/')]):
				os.mkdir(temp[:temp.rfind('/')])
		
		# Read the profiles and models
		print("-------> Read and Write Profiles ...")
		stokes = p.profile_stk(0,0,0)
		stokes.read_results_MC(path, tasks, "best.per")
		stokes.write(f"{os.path.join(path,conf['inv_out'])}{d.end_stokes}")

		if conf['chi2'] == "":
			del stokes

		print("-------> Read and Write Models ...")
		models = m.model_atm(conf["num"],1,0)	# Model
		errors = m.model_atm(conf["num"],1,0)	# Error
		guess  = m.model_atm(conf["num"],1,0)   # Best guess model
		
		models.read_results(tasks, "best.mod", path, int(conf['num']), 1)
		errors.read_results(tasks, "best.err", path, int(conf['num']), 1)
		guess.read_results(tasks, d.best_guess, path, int(conf['num']), 1)
		
		models.write(f"{os.path.join(path,conf['inv_out']+d.end_models)}")
		errors.write(f"{os.path.join(path,conf['inv_out']+d.end_errors)}")
		guess. write(f"{os.path.join(path,conf['inv_out']+d.best_guess_file)}")

		del models
		del errors
		del guess

		if conf['chi2'] == "1":
			print("-------> Compute χ² (this might take a while) ...")
			chi2 = c.chi2_stk(stokes.nx,stokes.ny)

			# Number of Nodes in the last step
			num_of_nodes = 0
			if conf['nodes_temp'] != "":
				num_of_nodes+= int(conf['nodes_temp'].split(",")[-1])
			if conf['nodes_magn'] != "":
				num_of_nodes+= int(conf['nodes_magn'].split(",")[-1])
			if conf['nodes_vlos'] != "":
				num_of_nodes+= int(conf['nodes_vlos'].split(",")[-1])
			if conf['nodes_gamma'] != "":
				num_of_nodes+= int(conf['nodes_gamma'].split(",")[-1])
			if conf['nodes_phi'] != "":
				num_of_nodes+= int(conf['nodes_phi'].split(",")[-1])
			
			# Compute chi2
			chi2.compute(obs, stokes, [float(i) for i in conf["weights"]], num_of_nodes)
			print("-------> Total χ² = %.3f" % chi2.total)
			chi2.write(os.path.join(path, conf['inv_out'] + d.end_chi2))

			del obs
			del num_of_nodes
			del stokes			
			del chi2


		if not debug:
			print("-------> Delete created SIR files")
			for i in range(conf['num']):
				shutil.rmtree(os.path.join(path,tasks['folders'][i]))
			
			os.remove(os.path.join(path,d.inv_trol_file))
			os.remove(os.path.join(path,d.Grid))
		else:
			print("-------> No created files are deleted. Debug option active.")

		# Print needed time
		end = time.time()
		Time = str(datetime.timedelta(seconds=end-start)).split(".")[0]
		print(f"-------> Finished in {Time}.")

	comm.barrier()

	return

"""
*****************************************************************************
*							TWO COMPONENT INVERSION							*
*									START									*
*****************************************************************************
"""

def inversion_2c(conf : dict, comm, rank : int, size : int, MPI, debug=False,progress=True):
	"""
	Performs the inversion of all the models for the 2 component inversion.

	Parameters
	----------
	conf : dict
		Dictionary with all information from the config file
	comm : Intercomm
		MPI communicator
	rank : int
		Number of actual process
	size : int
		Number of available processes
	MPI : library
		MPI library imported as 'from mpi4py import MPI'
	debug : bool,optional
		Execute in debugging mode
	progress : bool,optional
		Print a progress bar
	
	Returns
	-------
	None

	"""
	start = time.time()
	###################################################
	#		READ PARAMETERS FROM CONFIG FILE		#	
	###################################################
	# Define parameters for easier understanding
	path = conf["path"]
	model1 = conf["model1"]
	model2 = conf["model2"]

	line_file = conf['line'] # line file

	# Read inversion stuff
	Map = conf['map']
	if rank == 0:
		print(f"-------> Selected Map ({Map[0]},{Map[1]},{Map[2]},{Map[3]})")
	

	abundance_file = conf['abundance'] # Abundance file	
	end='.per' # ending of the file for the profile

	# Write psf function, if needed
	if rank == 0:
		if conf['psf'] != '':
			#print("-------> Spectral PSF is used")
			if conf['psf'].split(" ")[0] == "gauss":
				print(r"-------> Create Gaussian spectral PSF with σ = " + conf['psf'].replace("gauss ","") + " mÅ")
				sir.write_gauss_psf(float(conf['psf'].replace("gauss ","")), os.path.join(path,d.psf))
			else:
				print(r"-------> Use spectral PSF file {conf['psf']}")
				if(conf['psf'] != d.psf):
					shutil.copy(os.path.join(conf['path'],conf['psf']),os.path.join(conf['path'],d.psf))

	# Create guess from bin file if wanted
	if conf["guess1"] != '':
		guess1 = None
		if rank == 0:
			guess1 = m.read_model(os.path.join(path, conf["guess1"].split(" ")[0]))
			print(f"-------> File {conf['guess1'].split(' ')[0]} used as initial guess/base model 1")
			if guess1.nx > Map[1]-Map[0]+1 or guess1.ny > Map[3]-Map[2]+1:
					print("[Warn]   The shapes of the initial guess/base model 1 are bigger than the map in the config.")
			if guess1.nx < Map[1]-Map[0]+1 or guess1.ny < Map[3]-Map[2]+1:
					print(f"[ERROR]  The shapes ({guess1.nx},{guess1.ny}) of the initial guess/base model 1 are not compatible with the map ({Map[1]-Map[0]+1},{Map[3]-Map[2]+1}) in the config.")
					return
			
			# Small variation around the guess
			if(len(conf["guess1"].split(" ")) > 1):
				factor1 = float(conf["guess1"].split(" ")[1])
				print(f"-------> Create small variation around the initial guess with the factor {factor1}")
				guess1 = g.create_small_guess(guess1, conf["random_pars"], factor1)
		
		guess1 = comm.bcast(guess1, root=0)

	if conf["guess2"] != '':
		guess2 = None
		if rank == 0:
			guess2 = m.read_model(os.path.join(path, conf["guess2"].split(" ")[0]))
			print(f"-------> File {conf['guess2'].split(' ')[0]} used as initial guess/base model 2")
			if guess2.nx > Map[1]-Map[0]+1 or guess2.ny > Map[3]-Map[2]+1:
					print("[Warn]   The shapes of the initial guess/base model 2 are bigger than the map in the config.")
			if guess2.nx < Map[1]-Map[0]+1 or guess2.ny < Map[3]-Map[2]+1:
					print(f"[ERROR]  The shapes ({guess2.nx},{guess2.ny}) of the initial guess/base model 2 are not compatible with the map ({Map[1]-Map[0]+1},{Map[3]-Map[2]+1}) in the config.")
					return
			
			# Small variation around the guess
			if(len(conf["guess2"].split(" ")) > 1):
				factor2 = float(conf["guess2"].split(" ")[1])
				print(f"-------> Create small variation around the initial guess with the factor {factor2}")
				guess2 = g.create_small_guess(guess2, conf["random_pars"], factor2)
		
		guess2 = comm.bcast(guess2, root=0)

	

	# Print out randomisation setting
	if rank == 0:
		# Write which parameters are randomised
		if conf["random_guess"] > 0 or (conf["guess1"] != '' and len(conf["guess1"].split(" ")) > 1) or (conf["guess2"] != '' and len(conf["guess2"].split(" ")) > 1):
			random_pars = conf["random_pars"]
			print("-------> Parameters to be randomised: ", end = '')
			print_out = ''
			if 'T' in random_pars:
				print_out += "Temperature, "
			if 'Pe' in random_pars:
				print_out += "El. Pressure, "
			if 'vmicro' in random_pars:
				print_out += "Microturbulent Velocity, "
			if 'B' in random_pars:
				print_out += "Magnetic Field, "
			if 'vlos' in random_pars:
				print_out += "Line-of-Sight Velocity, "
			if 'gamma' in random_pars:
				print_out += "Inclination, "
			if 'phi' in random_pars:
				print_out += "Azimuth, "
			if 'z' in random_pars:
				print_out += "Height, "
			if 'Pg' in random_pars:
				print_out += "Gas Pressure, "
			if 'rho' in random_pars:
				print_out += "Density, "
			print(print_out[0:-2] + ".")
			if conf['random_guess'] == 1:
				print(f"-------> Create {conf['random_guess']} random guess per pixel.")
			else:
				print(f"-------> Create {conf['random_guess']} random guesses per pixel.")
		else:
			if conf["guess1"] == '':
				print(f"-------> Use Model '{conf['model1']}' as initial guess 1.")
			if conf["guess2"] == '':
				print(f"-------> Use Model '{conf['model2']}' as initial guess 2.")

	if rank == 0:
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "_*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)

	if rank == 0:
		fills = [float(i) for i in conf["fill"].split(",")]
		if (fills[0] + fills[1] != 1.0):
			print(f"-------> WARNING: Filling factors do not add up to 1.0 ({fills[0]} + {fills[1]} = {'%.2f' % (fills[0] + fills[1])})")

	# Load and scatter data => Saving memory and time
	if conf['chi2'] != "":
		stk, tasks, obs = scatter_data(conf, comm, rank, size, True)
	else:
		stk, tasks = scatter_data(conf, comm, rank, size, False)

	# Write the control file with the information from the config file
	if rank == 0:
		print(f"-------> Write control and grid file")
		sir._write_control_2c(os.path.join(path,d.inv_trol_file), conf)
		sir.write_grid(conf, os.path.join(path,d.Grid), stk.wave)


	if rank == 0:
		print("[STATUS] Start Computing Inversions ...")

	performed_models = 0  # Counts how many models are performed in each rank
	total_jobs = 1  # Total performed jobs across all processes
	temp_tasks = sir.create_task_folder_list(Map) # Structure tasks
	max_jobs = len(temp_tasks['folders'])  # For comm.allreduce function
	del temp_tasks

	

	# Root process initializes the progress bar
	if rank == 0 and progress:
		pbar = tqdm(total=max_jobs, desc="Overall Progress", file=sys.stdout)

	comm.barrier()

	for i in range(0, len(tasks['folders'])):
		#########################################
		#	Create the folder and copy stuff	#
		#########################################
		# Read x and y position, create task folder
		task_folder = os.path.join(path, tasks['folders'][i])
		
		x = tasks['x'][i]
		y = tasks['y'][i]
		os.makedirs(task_folder, exist_ok=True)
					
		# Write the data from the cube into a profile file for SIR
		stk.write_profile(os.path.join(path, task_folder, d.profile_obs), i, 0, os.path.join(conf['path'],d.Grid))
		
		# Copy stuff for the inversion
		if conf['psf'] != '':
			sir_files_src = [d.inv_trol_file, model1,       model2, "sir.x", line_file, d.Grid, abundance_file, d.psf]
			sir_files_dst = [d.inv_trol_file, model1 + "1", model2 + "2", "sir.x", line_file, d.Grid, abundance_file, d.psf]
		else:
			sir_files_src = [d.inv_trol_file, model1,       model2, "sir.x", line_file, d.Grid, abundance_file]
			sir_files_dst = [d.inv_trol_file, model1 + "1", model2 + "2", "sir.x", line_file, d.Grid, abundance_file]

		for i in range(len(sir_files_src)):
			shutil.copy(os.path.join(path, sir_files_src[i]), os.path.join(task_folder, sir_files_dst[i]))
		
		
		# Use a guess file for model 1
		if conf["guess1"] != '':
			# Write the new initial model from here:
			x1 = x - Map[0] # x position taking the Map into account
			y1 = y - Map[2] # y position taking the Map into account
			guess1.write_model(os.path.join(task_folder,model1 + "1"),x1,y1)
			
		# Use a guess file for model 1
		if conf["guess2"] != '':
			# Write the new initial model from here:
			x1 = x - Map[0] # x position taking the Map into account
			y1 = y - Map[2] # y position taking the Map into account
			guess1.write_model(os.path.join(task_folder,model2 + "2"),x1,y1)

		# Create random guesses and select best value and perform inversion
		execute_inversion_2c(conf, task_folder, rank)

		performed_models += 1
	
		# Do not do allreduce for the last step as the code does not move on from here
		if total_jobs < (max_jobs - max_jobs % size):
			# Root process updates the progress bar
			total_jobs = comm.allreduce(performed_models, op=MPI.SUM)
		
		# Update progres bar
		if rank == 0 and progress:
			pbar.n = total_jobs
			pbar.refresh()
		elif rank == 0:
			print(f"\rTotal Jobs: {total_jobs}/{max_jobs}", end="")

	comm.barrier()
		
	if rank == 0 and progress:
		pbar.n = max_jobs
		pbar.refresh()
		pbar.close()
	elif rank == 0:
			print(f"\rTotal Jobs: {total_jobs}/{max_jobs}")

	########################################################
	#	Read all the results and put it into npy files	#
	########################################################
	if rank == 0:
		
		print("[STATUS] Gathering results...")

		# Redefine tasks as now all the tasks are read
		tasks = sir.create_task_folder_list(Map) # Structure tasks

		#  Create directory if / exists in inv_out in config
		if '/' in conf['inv_out']:
			temp = os.path.join(path, conf["inv_out"])
			if not exists(temp[:temp.rfind('/')]):
				os.mkdir(temp[:temp.rfind('/')])

		print("-------> Read and Write Profiles ...")
		stokes_inv = p.profile_stk(stk.nx, stk.ny, stk.nw)

		stokes_inv.wave = stk.wave # Copy wavelength positions
		stokes_inv.read_results(tasks, f"best.per", path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)

		stokes_inv.write(os.path.join(path,conf['inv_out']) + d.end_stokes)

		if conf["chi2"] == "":
			del stokes_inv

		print("-------> Read and Write Models ...")
		# Create shapes of the arrays which are filled and saved later
		log_tau, _,_,_,_,_,_,_,_,_,_ = sir.read_model(f"{os.path.join(path,tasks['folders'][0])}/best1.mod")

		models_inv1		= m.model_atm(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, nval=len(log_tau))
		models_inv2		= m.model_atm(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, nval=len(log_tau))
		errors_inv1		= m.model_atm(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, nval=len(log_tau))
		errors_inv2		= m.model_atm(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, nval=len(log_tau))
		best_guesses1	= m.model_atm(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, nval=len(log_tau))
		best_guesses2	= m.model_atm(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, nval=len(log_tau))

		del log_tau

		models_inv1.read_results(tasks, 'best1.mod', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		models_inv2.read_results(tasks, 'best2.mod', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		errors_inv1.read_results(tasks, 'best1.err', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		errors_inv2.read_results(tasks, 'best2.err', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		best_guesses1.read_results(tasks, d.best_guess1, path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		best_guesses2.read_results(tasks, d.best_guess2, path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		
		models_inv1.write(os.path.join(path,conf['inv_out'] + d.end_models1))
		models_inv2.write(os.path.join(path,conf['inv_out'] + d.end_models2))
		errors_inv1.write(os.path.join(path,conf['inv_out'] + d.end_errors1))
		errors_inv2.write(os.path.join(path,conf['inv_out'] + d.end_errors2))
		best_guesses1.write(os.path.join(path,conf['inv_out'] + d.best_guess1_file))
		best_guesses2.write(os.path.join(path,conf['inv_out'] + d.best_guess2_file))

		del models_inv1
		del models_inv2
		del errors_inv1
		del errors_inv2
		del best_guesses1
		del best_guesses2

		if conf['chi2'] == "1":
			print("-------> Compute χ² (this might take a while) ...")
			chi2 = c.chi2_stk(stokes_inv.nx,stokes_inv.ny)

			# Number of Nodes in the last step
			num_of_nodes = 0
			if conf['nodes_temp1'] != "":
				num_of_nodes+= int(conf['nodes_temp1'].split(",")[-1])
			if conf['nodes_magn1'] != "":
				num_of_nodes+= int(conf['nodes_magn1'].split(",")[-1])
			if conf['nodes_vlos1'] != "":
				num_of_nodes+= int(conf['nodes_vlos1'].split(",")[-1])
			if conf['nodes_gamma1'] != "":
				num_of_nodes+= int(conf['nodes_gamma1'].split(",")[-1])
			if conf['nodes_phi1'] != "":
				num_of_nodes+= int(conf['nodes_phi1'].split(",")[-1])
			if conf["invert_fill"] != "1":
				num_of_nodes +=1
			
			# Compute chi2
			chi2.compute(obs, stokes_inv, [float(i) for i in conf["weights"]], num_of_nodes)
			print("-------> Total χ² = %.3f" % chi2.total)
			print("-------> Write χ²")
			
			chi2.write(os.path.join(path, conf['inv_out'] + d.end_chi2))
			del chi2
			del obs
			del num_of_nodes
			del stokes_inv

		
		if not debug:
			# Delete the folder
			print("-------> Delete created SIR files")
			for i in range(len(tasks['x'])):
				shutil.rmtree(os.path.join(path,tasks['folders'][i]))
			if conf['psf'] != '' and conf['psf'] != d.psf:
				os.remove(os.path.join(path,d.psf))
			os.remove(os.path.join(path,d.inv_trol_file))
			os.remove(os.path.join(path,d.Grid))
		else:
			print("-------> No created files are deleted. Debug option active.")

		# Print needed time
		end = time.time()
		Time = str(datetime.timedelta(seconds=end-start)).split(".")[0]
		print(f"-------> Finished in {Time}")

	comm.barrier()

	return
