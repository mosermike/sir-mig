#!/usr/bin/env python
"""
Perform inversion of many profiles with simple parallelization using mpi4py
"""

import os
import shutil
import numpy as np
import time
import glob
import sys
import sir
from create_random_guess import create_guesses_1c as create_guesses
from os.path import exists
import datetime
import definitions as d
import profile_stk as p
import model as m

def help():
	"""
	Help page
	"""	
	print("inversion - Executes the inversion")
	print("Usage: python inversion.py [OPTION]")
	print()
	print("1:  Config file")
	sys.exit()
	
################################################################################

def pprint(text, rank):
	"""
	Print out status information for rank 0

	Parameter
	---------
	text : str
		String which is printed
	rank : integer
		The rank which executes this function

	Return
	------
	None
	"""
	if rank == 0:
		print(text)

################################################################################

def x_y_add_zeros(x,y):
	"""
	Adds zeros so that the returning strings have 4 letters

	Parameter
	---------
	x : float
		x position
	y : float
		y position

	Return
	------
	x as a string of 4 letters
	y as a string of 4 letters
	"""
	if x < 10:
		x_str = "000" + str(x)
	elif x < 100:
		x_str = "00" + str(x)
	elif x < 1000:
		x_str = "0" + str(x)
	else:
		x_str = str(x)
	if y < 10:
		y_str = "000" + str(y)
	elif y < 100:
		y_str = "00" + str(y)
	elif y < 1000:
		y_str = "0" + str(y)
	else:
		y_str = str(y)

	return x_str, y_str


def create_task_folder_list(num):
	"""
	Creates a list which folders should be created and executed. This is done so
	that the inversion itself can be executed linearly to make use of all cores.

	Parameter
	---------
	num : int
		Number of models

	Return
	------
	Dictionary with all the names of the task folders, x and y position
	"""
	# Create arrays
	tasks = []
	xs = []
	ys = []
	
	y = 0
	# Determine task folder names
	for x in range(num):
		x_str, y_str = x_y_add_zeros(x,0)

		tasks.append(d.task_start + x_str + "_" + y_str)
		xs.append(x)
		ys.append(y)
	Dict = {
			'folders' : tasks,
			'x' : np.array(xs),
			'y' : np.array(ys),
	
	}

	return Dict
	
def scatter_data(conf, comm, rank, size):
	"""
	Scatters the data equally into all the ranks

	Parameter
	---------
	conf : Dict
		Config information
	path : string
		Path where we are
	"""
	path = conf["path"]
	

	if rank == 0:
		print("[STATUS] Load and scatter data ...")
		
		if "--no-noise" in sys.argv:
			if rank == 0:
				print("-------> No noise flag used")
				print("-------> Use synthesis profiles")
			stk = p.read_profile(os.path.join(path,conf["syn_out"]))
		else:
			stk = p.read_profile(os.path.join(path,conf["noise_out"]))

		tasks = create_task_folder_list(conf["num"])

		# Create one data cube
		stki = stk.stki
		stkq = stk.stkq
		stku = stk.stku
		stkv = stk.stkv
		lines = stk.indx
		waves = stk.wave

		del stk # Free Memory

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
			stki_chunks[-1] = stki[-remainder-chunk_size:]
			stkq_chunks[-1] = stkq[-remainder-chunk_size:]
			stku_chunks[-1] = stku[-remainder-chunk_size:]
			stkv_chunks[-1] = stkv[-remainder-chunk_size:]
			tasks_chunks[-1] = tasks['folders'][-remainder-chunk_size:]
			num_chunks[-1] = tasks['x'][-remainder-chunk_size:]
	else:
		stki_chunks = None
		stkq_chunks = None
		stku_chunks = None
		stkv_chunks = None
		num_chunks = None
		tasks_chunks = None
		waves = None
		lines = None

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
	stk = p.Profile(stki_chunk.shape[0],stki_chunk.shape[1],stki_chunk.shape[2])
	stk.data_cut = True # Already cut before
	stk.indx = lines
	stk.stki[:,:] = stki_chunk
	stk.stkq[:,:] = stkq_chunk
	stk.stku[:,:] = stku_chunk
	stk.stkv[:,:] = stkv_chunk
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

	return stk, tasks

################################################################################

def execute_inversion(conf, task_folder):
	"""
	Executes inversion and creates if needed random guesses. Make sure that
	when this function is entered, the os is in the right directory!
	
	Parameter
	---------
	config : dict
		Dictionary containing the configuration of the simulation

	Return
	------
	None
	"""
	# Define parameters for simplicity reasons
	shutil.copy(conf['model'], d.model_inv)
	model = d.model_inv
	guess1 = d.model_inv.replace(".mod", "")  # For simplicity
	cycles = conf["cycles"]
	chi_file = d.inv_trol_file[:d.inv_trol_file.rfind('.')] + ".chi"
	
	# Remove chi2 result file in case it exists
	if exists(chi_file):
		os.remove(chi_file)

	###################
	# START INVERSION #
	###################
	# Determine whether random guesses are wished or not
	if conf["random_guess"] > 0:
		# Initialize the chi2 map
		chi2 = np.zeros(conf["random_guess"])

		# Perform inversion for each guess model and copy files

		for i in range(conf["random_guess"]):
			do_inversion = True
			it = 0
			while do_inversion and it < 50:
				# Create New Guess
				create_guesses(conf, output="./", number=i+1)
				# Copy to the model
				shutil.copy(d.model + str(i+1) + ".mod", d.model_inv)
				# Execute inversion
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")
				# Sometimes it does not converge and when only 8 cols are used, the programme could crash
				# => results into NaNs
				# Other bugs could result it to not finish properly => Check dimension of inv.chi file
				# Create a new random guess and rerun
				if exists(chi_file):  # Chi file exists
					if np.genfromtxt('inv.chi').shape == (int(cycles), 2) or int(cycles) == 1:  # All the cycles finished
						break
					else:
						os.remove(chi_file)
				else: # No inv.chi file generated => There might be a problem with sir.x => Break loop after it is greater than 50
					it += 1

			# If it is greater than 50, there might be something wrong with sir.x
			if it >= 50:
				print("[ERROR] Check your sir.x file and the log file in the .task folders. There might be a problem with sir.x")
				# Print last log entry
				with open('inv.log') as f:
					for line in f:
						pass
					print("[LAST LOG ENTRY]: ", line)
				sys.exit()
					
			chi = sir.read_chi2(chi_file)

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
			
			# Remove not needed files
			for j in range(int(cycles)-1):
				os.remove(f"{guess1}_{j+1}.mod")
				os.remove(f"{guess1}_{j+1}.err")
				os.remove(f"{guess1}_{j+1}.per")

		# All inversions failed => Redo inversion until one does not fail
		if chi2.min() > 1e8:
			it = 0  # Number of iterations until sir converges
			chi2_best = 0.0
			while chi2_best == 0.0:
				# Create New Guess
				create_guesses(conf, output="./", number=i+1)

				# Copy to the model
				shutil.copy(guess1 + str(i+1) + ".mod", model)

				# Execute inversion again
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Read chi2 value
				chi2_best = sir.read_chi2(chi_file)
	
				# +1 iteration
				it += 1
			# Remove not needed files
			for i in range(int(cycles)-1):
				os.remove(f"{guess1}_{i+1}.mod")
				os.remove(f"{guess1}_{i+1}.err")
				os.remove(f"{guess1}_{i+1}.per")

			shutil.copy(f"{guess1}.mod", f"{d.best_guess}")  # Copy best guess model

			# Warn if the repetition of the inverse needed to be done more than 20 times
			if it > 10:
				print(f"\nNote that the while loop run in {task_folder[task_folder.rfind('/')+1:]} took {it} iterations because chi2 was 0.0...")

		# Determine best fit and copy stuff
		else:

			chi2_best = chi2.argmin()
			shutil.copy(f"inv.chi__{str(chi2_best+1)}", chi_file)
			shutil.copy(f"{guess1}__{str(chi2_best+1)}.mod", f"{guess1}_{cycles}.mod")
			shutil.copy(f"{guess1}__{str(chi2_best+1)}.err", f"{guess1}_{cycles}.err")
			shutil.copy(f"{guess1}__{str(chi2_best+1)}.per", f"{guess1}_{cycles}.per")
			shutil.copy(f"{d.model}{str(chi2_best+1)}.mod", f"{d.best_guess}")  # Copy best guess model

	else:
		# Perform inversion once and read chi2 value
		os.system(f"echo {d.inv_trol_file} | ./sir.x > /dev/null")
		shutil.move(model, d.best_guess)
		# Remove not needed files
		for i in range(int(cycles)-1):
			os.remove(f"{guess1}_{i+1}.mod")
			os.remove(f"{guess1}_{i+1}.err")
			os.remove(f"{guess1}_{i+1}.per")

	return


def inversion(conf, comm, rank, size, MPI):
	"""
	Performs the inversion of all the models.

	Parameter
	---------
	conf : dict
		Dictionary with all information from the config file
	comm, rank, size, MPI : MPI variables
		Variables defined as
		 - comm = MPI.COMM_WORLD
		 - rank = comm.Get_rank()
		 - size = comm.Get_size()
		 - MPI = MPI library imported as from mpi4py import MPI
	

	Return
	------
	None

	"""
	####################################
	# READ PARAMETERS FROM CONFIG FILE #
	####################################
	# Define parameters for easier understanding
	path = conf["path"]
	model = conf["model"]
	line_file = conf['line']  # line file
	
	# Write the control file with the information from the config file
	if rank == 0:
		sir.write_control_mc(os.path.join(conf['path'], d.inv_trol_file), conf)
	abundance_file = conf['abundance']  # Abundance file
	

	# Create guess from npy file if wanted
	if conf["guess"] != '':
		guess = np.load(os.path.join(path, conf["guess"]))
		pprint(f"-------> File {conf['guess']} used as initial guess/base model", rank)

	####################
	# CREATE GRID FILE #
	####################
	if rank == 0:
		# Write Grid file based on the chosen wavelength ranges in the config file
		sir.write_grid_mc(conf, os.path.join(path, d.Grid))
		# Write which parameters are randomised
		if conf["random_guess"] > 0:
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
		else:
			pprint(f"-------> Use Base Model '{conf['model']}' as initial guess.", rank=rank)

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
	pprint("[STATUS] Start Inversions", rank=rank)
	performed_models = 0  # Counts how many models are performed
	total_jobs = 1  # Total performed jobs across all processes
	max_jobs = conf['num']  # For comm.allreduce function


	stk, tasks = scatter_data(conf, comm, rank, size)

	start_time = time.time()
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
		
		# Create guess from npy file if wanted
		# Copy the data from the imported numpy array to the model in the task folder
		if conf["guess"] != '':
			# Write the new initial model from here:
			sir.write_model(os.path.join(task_folder, model), d.header,
							guess[i, 0], guess[i, 1], guess[i, 2], guess[i, 3],
							guess[i, 4], guess[i, 5], guess[i, 6], guess[i, 7],
							guess[i, 8], guess[i, 9], guess[i, 10]
						)

		###############################
		# 	Perform inversion		#
		###############################
		# Create random guesses and select best value and perform inversion
		os.chdir(task_folder)
		
		execute_inversion(conf, task_folder)
		
		performed_models += 1
	
		# Do not do allreduce for the last step as the code does not move on from here
		if total_jobs < (max_jobs - max_jobs % size):
			total_jobs = comm.allreduce(performed_models, op=MPI.SUM)

		# Print the total number of finished jobs on the root process
		if rank == 0:
			elapsed_time = time.time() - start_time
			remaining_time = elapsed_time * max_jobs/total_jobs - elapsed_time
			remaining_time = str(datetime.timedelta(seconds=remaining_time)).split(".")[0]  # Convert time into clock format
			print(f"\rFinished Jobs: [{total_jobs}/{max_jobs}], Remaining time {remaining_time}", end='', flush=False)

		os.chdir("../")  # In case a relative path is used in the config

	comm.barrier()
		
	####################################
	# Print out time and finished jobs #
	####################################
	if rank == 0:
		remaining_time = str(datetime.timedelta(seconds=0)).split(".")[0]
		total_time = time.time() - start_time  # Total time
		total_time = str(datetime.timedelta(seconds=total_time)).split(".")[0]  # Remove ms
		print(f"\r-------> Finished {max_jobs} Jobs in {total_time}.                          ")

	##################################################
	# Read all the results and put it into npy files #
	##################################################
	os.chdir(path)
	if rank == 0:
		start = time.time()

		tasks = create_task_folder_list(conf["num"])

		# Read the profiles and models
		stokes = p.Profile(0,0,0)
		stokes.read_results_MC(path, tasks, f"{d.model_inv.replace('.mod','')}_{conf['cycles']}.per")
		
		models = m.Model(conf["num"],1,0)	# Model
		errors = m.Model(conf["num"],1,0)	# Error
		guess  = m.Model(conf["num"],1,0)   # Best guess model
		
		models.read_results(tasks, f"{d.model_inv.replace('.mod','')}_{conf['cycles']}.mod", path, int(conf['num']), 1)
		errors.read_results(tasks, f"{d.model_inv.replace('.mod','')}_{conf['cycles']}.err", path, int(conf['num']), 1)
		guess.read_results(tasks, d.best_guess, path, int(conf['num']), 1)
		
		chi2 = sir.read_chi2s(conf, tasks)
		
		stokes.write(f"{os.path.join(path,conf['inv_out'])}{d.end_stokes}")
		
		models.write(f"{os.path.join(path,conf['inv_out'])}{d.end_models}")
		errors.write(f"{os.path.join(path,conf['inv_out'])}{d.end_errors}")
		guess.write(f"{os.path.join(path,d.best_guess.replace('.mod','.bin'))}")
		np.save(f"{conf['chi2']}", chi2)

		for i in range(conf['num']):
			shutil.rmtree(os.path.join(path,tasks['folders'][i]))

		# Print needed time
		end = time.time()
		Time = str(datetime.timedelta(seconds=end-start)).split(".")[0]
		print(f"-------> Finished in {Time}.")

	comm.barrier()

	return
