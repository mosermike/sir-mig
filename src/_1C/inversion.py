#!/usr/bin/env python
"""
Perform inversion of many profiles with simple parallelization using mpi4py
"""

import os
import shutil
import numpy as np
import time
import sys
import sir
import obs
import glob
import model as m
import profile_stk as p
from _1C.create_random_guess import create_guesses
from os.path import exists
import datetime
import definitions as d

def x_y_add_zeros(x, y):
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


def create_task_folder_list(Map):
	"""
	Creates a list which folders should be created and executed. This is done so
	that the inversion itself can be executed linearly to make use of all cores.

	Parameter
	---------
	Map : numpy array
		1x4 array containing the limits in x and y of the data

	Return
	------
	Dictionary with all the names of the task folders, x and y position
	"""
	# Create arrays
	tasks = []
	xs = []
	ys = []
	
	# Determine task folder names
	for x in range(Map[0], Map[1]+1):
		for y in range(Map[2], Map[3]+1):
			x_str, y_str = x_y_add_zeros(x, y)

			tasks.append(d.task_start + x_str + "_" + y_str)
			xs.append(x)
			ys.append(y)
	Dict = {
			'folders': tasks,
			'x': np.array(xs),
			'y': np.array(ys),
	
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
	tasks = create_task_folder_list(conf["map"])
	if rank == 0:
		print("[Status] Load and scatter data ...")
		stk = p.read_profile(os.path.join(path,conf["cube_inv"]))
		stk.cut_to_wave(conf["range_wave"]) # Cut wavelength file to the wished area
		# Create one data cube
		stki = stk.stki.reshape(-1, stk.nw) # Flatten Stokes I corresponding to the tasks list
		stkq = stk.stkq.reshape(-1, stk.nw) # Flatten Stokes Q corresponding to the tasks list
		stku = stk.stku.reshape(-1, stk.nw) # Flatten Stokes U corresponding to the tasks list
		stkv = stk.stkv.reshape(-1, stk.nw) # Flatten Stokes V corresponding to the tasks list		
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
		x_chunks = [tasks['x'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]
		y_chunks = [tasks['y'][i*chunk_size:(i+1)*chunk_size] for i in range(size)]

		# If there is a remainder, distribute the remaining data to the last process
		if remainder:
			stki_chunks[-1] = stki[-remainder-chunk_size:]
			stkq_chunks[-1] = stkq[-remainder-chunk_size:]
			stku_chunks[-1] = stku[-remainder-chunk_size:]
			stkv_chunks[-1] = stkv[-remainder-chunk_size:]
			tasks_chunks[-1] = tasks['folders'][-remainder-chunk_size:]
			x_chunks[-1] = tasks['x'][-remainder-chunk_size:]
			y_chunks[-1] = tasks['y'][-remainder-chunk_size:]
	else:
		stki_chunks = None
		stkq_chunks = None
		stku_chunks = None
		stkv_chunks = None
		x_chunks = None
		y_chunks = None
		tasks_chunks = None
		waves = None

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
	stk = p.Profile(stki_chunk.shape[0],1,nw=stki_chunk.shape[1])
	stk.data_cut = True # Already cut before
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
	chi2 : float
		Best chi2 value of the inversion(s)
	"""
	# Define parameters for simplicity reasons
	model = conf['model']
	guess1 = model.replace(".mod", "")  # For simplicity
	cycles = conf["cycles"]
	chi_file = d.inv_trol_file[:d.inv_trol_file.rfind('.')] + ".chi"

	#########################
	# START INVERSION #
	#########################
	# Determine whether random guesses are wished or not
	if conf["random_guess"] > 0:
		# Initialize the chi2 map
		chi2 = np.zeros(conf["random_guess"])

		# Create random guesses depending on the configurations
		create_guesses(conf, output="./")

		# Perform inversion for each guess model and copy files

		for i in range(conf["random_guess"]):
			it = 0
			while it < 50: # stop (repeating) inversion when chi2 file is correctly generated => inversion finished correctly or more than 50 iterations
				# Create New Guess
				create_guesses(conf, output="./", number=i+1)
				# Copy to the model
				shutil.copy(guess1 + str(i+1) + ".mod", model)
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
				print("[ERROR] Check your sir.x file and the log file in the .task folders. There might be a problem with sir.x")
				# Print last log entry
				with open('inv.log') as f:
					for line in f:
						pass
					print("[LAST LOG ENTRY]: ", line)
				sys.exit()
			# Read the chi2 file
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

			# Remove not needed files to make space for the next inversion
			for j in range(int(cycles)-1):
				os.remove(f"{guess1}_{j+1}.mod")
				os.remove(f"{guess1}_{j+1}.err")
				os.remove(f"{guess1}_{j+1}.per")

		# All inversions failed => Redo inversion until one does not fail
		if chi2.min() > 1e8:
			it = 0 # Number of iterations until sir converges
			chi2_best = 0.0
			while chi2_best == 0.0:
				# Copy the guess which did not work
				shutil.copy(model, "bad_guess" + str(i+1) + ".mod")

				# Create New Guess
				create_guesses(conf, output="./", number=i+1)

				# Copy to the model
				shutil.copy(guess1 + str(i+1) + ".mod",model)

				# Execute inversion again
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Read chi2 value
				chi2_best = sir.read_chi2(chi_file)
	
				# +1 iteration
				it += 1
			
			shutil.copy(f"{guess1}.mod", f"{d.best_guess}")  # Copy best guess model

			# Warn if the repetition of the inverse needed to be done more than 20 times
			if it > 10:
				print(f"\nNote that the while loop run in {task_folder[task_folder.rfind('/')+1:]} took {it} iterations because chi2 was 0.0...")

		# Determine best fit and copy stuff
		else:

			chi2_best = chi2.argmin()
			shutil.copy(f"inv.chi__{str(chi2_best+1)}",chi_file)
			shutil.copy(f"{guess1}__{str(chi2_best+1)}.mod", f"{guess1}_{cycles}.mod")
			shutil.copy(f"{guess1}__{str(chi2_best+1)}.err", f"{guess1}_{cycles}.err")
			shutil.copy(f"{guess1}__{str(chi2_best+1)}.per", f"{guess1}_{cycles}.per")
			shutil.copy(f"{guess1}{str(chi2_best+1)}.mod", f"{d.best_guess}") # Copy best guess model
			chi2_best = chi2[chi2_best]

	else:
		# Check if vmacro is different in the header of the model and if yes change it
		header = np.loadtxt(conf['model'], max_rows=1)
		if header[0] != float(conf['vmacro']):
			header[0] = float(conf['vmacro'])
			tem = np.loadtxt(conf['model'], skiprows=1)
			if len(tem) > 8:
				sir.write_model(conf['model'], header, tem[0], tem[1], tem[2], tem[3], tem[4], tem[5], tem[6], tem[7], tem[8], tem[9], tem[10])
			else:
				sir.write_model(conf['model'], header, tem[0], tem[1], tem[2], tem[3], tem[4], tem[5], tem[6], tem[7])
		
		# Perform inversion once and read chi2 value
		os.system(f"echo {d.inv_trol_file} | ./sir.x > /dev/null")
		chi2_best = sir.read_chi2(chi_file)
		shutil.copy(model, d.best_guess)
		
		# Remove not needed files to make space for the next inversion
		for j in range(int(cycles)-1):
			os.remove(f"{guess1}_{j+1}.mod")
			os.remove(f"{guess1}_{j+1}.err")
			os.remove(f"{guess1}_{j+1}.per")

	return chi2_best


def inversion(conf, comm, rank, size, MPI):
	"""
	Performs the inversion of all the models.

	Parameter
	---------
	comm, rank, size : MPI variables
		Variables defined as
		 - comm = MPI.COMM_WORLD
		 - rank = comm.Get_rank()
		 - size = comm.Get_size()
	config : dict
		Dictionary with all information from the config file

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
	model1 = model.replace('.mod','') # For simplicity reasons defined
	
	#output =  os.path.join(path,conf["inv_out"]) # Output file
	cycles = conf['cycles'] # How many cycles
	line_file = conf['line'] # line file

	# Read inversion stuff
	Map = conf['map']
	
	# Write the control file with the information from the config file
	if rank == 0:
		print("-------> Write control and grid file")
		sir.write_control_1c(os.path.join(conf['path'],d.inv_trol_file), conf)
		# Write Grid file based on the chosen wavelength ranges in the config file
		stk = p.read_profile(os.path.join(path,conf["cube_inv"]))
		sir.write_grid(conf, stk.wave, os.path.join(path,d.Grid))
		del stk
	abundance_file = conf['abundance']  # Abundance file
	# Write psf function, if needed
	if rank == 1 or (size < 2 and rank == 0):
		if conf['psf'] != '':
			print("-------> Spectral PSF is used")
			if not exists(os.path.join(path, conf['psf'])):
				obs.write_psf(conf, os.path.join(path, conf['psf']))
	
	# Create guess from npy file if wanted
	if conf["guess"] != '':
		if rank == 0:
			print(f"-------> File {conf['guess']} used as initial guess/base model")
		guess = m.read_model(os.path.join(path,conf["guess"]))
		
	if rank == 2 or (size < 3 and rank == 0):
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)
	
	if rank == 0:
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
			print(f"-------> Create {conf['random_guess']} random guesses per pixel.")
		else:
			print(f"-------> Use Base Model '{conf['model']}' as initial guess.")

	########################
	# START INVERSION PART #
	########################
	performed_models = 0  # Counts how many models are performed
	total_jobs = 1  # Total performed jobs across all processes
	tasks = create_task_folder_list(Map) # Structure tasks
	max_jobs = len(tasks['folders'])  # For comm.allreduce function

	#########################
	# SCATTER AND LOAD DATA #
	#########################
	# Load and scatter data => Saving memory and time
	stk, tasks = scatter_data(conf, comm, rank, size)

	if rank == 0:
		print("[STATUS] Start Computing Inversions ...")
	start_time = time.time()
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
			sir_files = [d.inv_trol_file, model, "sir.x", line_file, d.Grid, abundance_file, conf['psf']]
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
			guess.write_model(os.path.join(task_folder, model), d.header, x1, y1)

		#####################
		# Perform inversion #
		#####################
		# Create random guesses and select best value and perform inversion
		os.chdir(task_folder)
		chi2_best = execute_inversion(conf, task_folder)
		
		#########################################
		# Check chi2 and print out informations #
		#########################################
		# If chi2 is not small, print out model number and do not delete the files
		#if d.chi2_verbose:
		#	if chi2_best > d.chi2_lim or chi2_best < 1e-2:
		#		chi2s = np.append(chi2s, chi2_best)
		#		chi2s_num = np.append(chi2s_num, f"{x}_{y}")

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

	wave = np.copy(stk.wave)
	del stk # Free Memory

	comm.barrier()
		
	####################################
	# Print out time and finished jobs #
	####################################
	if rank == 0:
		remaining_time = str(datetime.timedelta(seconds=0)).split(".")[0]
		print(f"\rFinished Jobs: [{max_jobs}/{max_jobs}], Remaining time {remaining_time}            ", end='', flush=False)
		total_time = time.time() - start_time  # Total time
		total_time = str(datetime.timedelta(seconds=total_time)).split(".")[0]  # Remove ms
		print(f"\n-------> Finished {max_jobs} Jobs in {total_time}.")

	##################################################
	# Read all the results and put it into npy files #
	##################################################
	os.chdir(path)
	if rank == 0:
		start = time.time()
		print("[STATUS] Gathering results...")

		# Redefine tasks as now all the tasks are read
		tasks = create_task_folder_list(Map) # Structure tasks

		# Create shapes of the arrays which are filled and saved later
		stokes_inv = p.Profile(0,0,0)
		stokes_inv.wave = wave # Copy wavelength positions
		models_inv = m.Model(0,0,0)
		errors_inv = m.Model(0,0,0)
		best_guesses = m.Model(0,0,0)

		print("-------> Read Models ...")
		models_inv.read_results(tasks, f'{model1}_{cycles}.mod', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		errors_inv.read_results(tasks, f'{model1}_{cycles}.err', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		best_guesses.read_results(tasks, d.best_guess, path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)

		# Correct for phi ambiguity (no difference between 0 and 180 measurable)
		#models_inv.correct_phi()

		print("-------> Read Profiles ...")
		stokes_inv.read_results(tasks, f"{model1}_{cycles}.per", path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)

		chi2 = np.zeros(shape=(Map[1]-Map[0]+1, Map[3]-Map[2]+1))

		# Collect data from task folders and delete the folder
		for i in range(len(tasks['x'])):
			# Get folder name and x,y positions
			folder = tasks['folders'][i]
			x = tasks['x'][i]
			y = tasks['y'][i]
			
			# Read chi2 file
			chi2[x-Map[0], y-Map[2]] = sir.read_chi2(f"{folder}/{d.inv_trol_file[:d.inv_trol_file.rfind('.')]}.chi")

			# Remove folder
			shutil.rmtree(folder)

		# Create directory if / exists in inv_out in config
		if '/' in conf['inv_out']:
			temp = os.path.join(path, conf["inv_out"])
			if not exists(temp[:temp.rfind('/')]):
				os.mkdir(temp[:temp.rfind('/')])
				
		# Save as npy files
		stokes_inv.write(os.path.join(path,conf['inv_out']) + d.end_stokes)
		models_inv.write(os.path.join(path,conf['inv_out']) + d.end_models)
		errors_inv.write(os.path.join(path,conf['inv_out']) + d.end_errors)
		best_guesses.write(os.path.join(path,d.best_guess.replace(".mod",".bin")))
		np.save(os.path.join(path,conf['chi2']),chi2)
		
		# Print needed time
		end = time.time()
		Time = str(datetime.timedelta(seconds=end-start)).split(".")[0]
		print(f"-------> Finished in {Time}")

	comm.barrier()

	return
