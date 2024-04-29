#!/usr/bin/env python
"""
Perform inversion of many profiles with simple parallelization using mpi4py
"""

import os, shutil, numpy as np
import time
import sys
import sir
import obs
from _2C.create_random_guess import create_guesses
from os.path import exists
import glob
import datetime
import definitions as d
import model as m
import profile_stk as p

def help():
	"""
	Help page
	"""	
	print("inversion - Executes the inversion")
	print("Usage: python inversion.py [OPTION]")
	print()
	print("1:  Config file")
	sys.exit()
	

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
	for x in range(Map[0],Map[1]+1):
		for y in range(Map[2], Map[3]+1):
			x_str, y_str = x_y_add_zeros(x,y)

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
	tasks = create_task_folder_list(conf["map"])
	if rank == 0:
		print("[Status] Load and scatter data ...")
		stk = p.read_profile(os.path.join(path,conf["cube_inv"]))
		stk.cut_to_wave(sir.angstrom_to_pixel(stk.wave, conf["range_wave"])) # Cut wavelength file to the wished area
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
	conf : Dictionary
		Dictionary containing the configuration of the simulation

	Return
	------
	chi2 : float
		Best chi2 value of the inversion(s)
	"""
	# Define parameters for simplicity reasons
	model1 = conf['model1']
	model2 = conf['model2']
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
				create_guesses(conf, output = "./", number = i+1)
				if not exists(f'{d.model1}' + str(i+1) + ".mod"):
					print(f'No model {i+1} in {task_folder} with the name {d.model1}{i+1}.mod')
					sys.exit()
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
				print("[ERROR] Check your sir.x file and the log file in the .task folders. There might be a problem with sir.x")
				# Print last log entry
				with open('inv.log') as f:
					for line in f:
						pass
					print("[LAST LOG ENTRY]: ", line)
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
				os.remove(f"{d.guess2.replace('.mod','')}_{i+1}.per")

		# All inversions failed => Redo inversion until one does not fail
		if chi2.min() > 1e8:
			it = 0 # Number of iterations until sir converges
			chi2_best = 0.0
			while chi2_best == 0.0:
				# Create New Guess
				create_guesses(conf, output = "./", number = 1)

				# Copy to the model
				shutil.copy("model1_1.mod",d.guess1)
				shutil.copy("model2_1.mod",d.guess2)

				# Execute inversion again
				os.system(f"echo {d.inv_trol_file} | ./sir.x >> inv.log 2>&1")

				# Read chi2 value
				chi2_best = sir.read_chi2(chi_file, task=task_folder)

				# +1 iteration
				it += 1
			
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
		# Perform inversion once and read chi2 value
		shutil.copy(model1,d.guess1)
		shutil.copy(model2,d.guess2)
		os.system(f"echo {d.inv_trol_file} | ./sir.x > /dev/null")
		chi2_best = sir.read_chi2(chi_file)
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
			os.remove(f"{d.guess2.replace('.mod','')}_{i+1}.per")

	return chi2_best
	
def inversion(conf, comm, rank, size, MPI):
	"""
	Performs the inversion of all the models.

	Parameter
	---------
	conf : dict
		Dictionary with all information from the config file
	comm, rank, size : MPI variables
		Variables defined as
		 - comm = MPI.COMM_WORLD
		 - rank = comm.Get_rank()
		 - size = comm.Get_size()
		 - MPI = MPI library imported as from mpi4py import MPI
	
	
	Return
	------
	None
	"""
	###################################################
	#		READ PARAMETERS FROM CONFIG FILE		#	
	###################################################
	# Define parameters for easier understanding
	path = conf["path"]
	model1 = conf["model1"]
	model2 = conf["model2"]

	output =  os.path.join(path,conf["inv_out"]) # Output file
	cycles = conf['cycles'] # How many cycles
	line_file = conf['line'] # line file

	# Read inversion stuff
	Map = conf['map']
	
	# Write the control file with the information from the config file
	if rank == 0:
		sir.write_control_2c(os.path.join(conf['path'],d.inv_trol_file), conf)
	abundance_file = conf['abundance'] # Abundance file	
	end='.per' # ending of the file for the profile

	# Write psf function, if needed

	if rank == 0:
		if conf['psf'] != '':
			print("-------> Spectral PSF is used")
			if not exists(os.path.join(path,conf['psf'])):
				obs.write_psf(conf, os.path.join(path,conf['psf']))

	# Create guess from npy file if wanted
	if conf["guess1"] != '':
			guess1 = m.read_model(os.path.join(path,conf["guess1"])) # Load data
			# Read Header Infos from the base model, if exists
			if exists(os.path.join(path,conf["model1"])):
				header1 = np.loadtxt(os.path.join(path,conf["model1"]), max_rows=1)
				header1 = f"    {header1[0]}    {header1[1]}    {header1[2]}"
			else:
				header1 = d.header_2c
			if rank == 0:
				print(f"-------> File {conf['guess1']} used as initial guess/base model 1")
				# Check if the shapes match
				if guess1.nx > Map[1]-Map[0]+1 or guess1.ny > Map[3]-Map[2]+1:
					print("[Warn]   The shapes of the initial guess/base model 1 are bigger than the map in the config.")
				elif guess1.nx < Map[1]-Map[0]+1 or guess1.ny < Map[3]-Map[2]+1:
					print(f"[ERROR] The shapes ({guess1.nx},{guess1.ny}) of the initial guess/base model 1 are not compatible with the map ({Map[1]-Map[0]+1},{Map[3]-Map[2]+1}) in the config.")
					sys.exit()
	# Create guess from npy file if wanted
	if conf["guess2"] != '':
			guess2 = m.read_model(os.path.join(path,conf["guess2"])) # Load data
			# Read Header Infos from the base model, if exists
			if exists(os.path.join(path,conf["model2"])):
				header2 = np.loadtxt(os.path.join(path,conf["model2"]), max_rows=1)
				header2 = f"    {header2[0]}    {header2[1]}    {header2[2]}"
			else:
				header2 = d.header_2c
		
			if rank == 0:
				print(f"-------> File {conf['guess2']} used as initial guess/base model 2")
				# Check if the shapes match
				if guess2.nx > Map[1]-Map[0]+1 or guess2.ny > Map[3]-Map[2]+1:
					print("[Warn]   The shapes of the initial guess/base model 2 are bigger than the map in the config.")
				if guess2.nx < Map[1]-Map[0]+1 or guess2.ny < Map[3]-Map[2]+1:
					print(f"[ERROR]  The shapes ({guess2.nx},{guess2.ny}) of the initial guess/base model 2 are not compatible with the map ({Map[1]-Map[0]+1},{Map[3]-Map[2]+1}) in the config.")
					return
	##########################
	#	CREATE GRID FILE	#
	##########################
	if rank == 0:
		# Write Grid file based on the chosen wavelength ranges in the config file
		sir.write_grid(conf, stk.wave, os.path.join(path,d.Grid))

	# Print out randomisation setting
	if rank == 0:
		# Write which parameters are randomised
		if conf["random_guess"] > 0:
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
		else:
			print(f"-------> Use Base Model '{conf['model']}' as initial guess.")

	if rank == 0:
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "_*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)

	comm.barrier() # Continues when all processes entered this function

	
	####################################
	#		START INVERSION PART	#
	####################################
	chi2s = np.array([]) # Save chi number for print out above chi_lim
	chi2s_num = np.array([], dtype=str) # Save task number for print out above chi_lim
	performed_models = 0 # Counts how many models are performed
	total_jobs = 1 # Total performed jobs across all processes
	tasks = create_task_folder_list(Map)
	max_jobs = len(tasks['folders']) # For comm.allreduce function

	# Load and scatter data => Saving memory and time
	stk, tasks = scatter_data(conf, comm)

	if rank == 0:
		print("[STATUS] Start Computing Inversions ...")
	start_time = time.time()
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
			sir_files = [d.inv_trol_file, model1, model2, "sir.x", line_file, d.Grid, abundance_file, conf['psf']]
		else:
			sir_files = [d.inv_trol_file, model1, model2, "sir.x", line_file, d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(task_folder, sir_file))
		
		# Create guess from npy file if wanted
		# Copy the data from the imported numpy array to the model in the task folder
		if conf["guess1"] != '':
			# Write the new initial model from here:
			x1 = x - Map[0] # x position taking the Map into account
			y1 = y - Map[2] # y position taking the Map into account
			guess1.write_model(os.path.join(task_folder,model1), header1,x1,y1)
			
		if conf["guess2"] != '':
			# Write the new initial model from here:
			x1 = x - Map[0] # x position taking the Map into account
			y1 = y - Map[2] # y position taking the Map into account
			guess1.write_model(os.path.join(task_folder,model2), header2,x1,y1)

		###############################
		# 	Perform inversion		#
		###############################
		# Create random guesses and select best value and perform inversion
		os.chdir(task_folder)
		chi2_best = execute_inversion(conf, task_folder)
		
		##############################################
		#	Check chi2 and print out informations	#
		##############################################
		# If chi2 is not small, print out model number and do not delete the files
		if d.chi2_verbose:
			if chi2_best > d.chi2_lim or chi2_best < 1e-2:
				chi2s = np.append(chi2s,chi2_best)
				chi2s_num = np.append(chi2s_num,f"{x}_{y}")
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

	comm.barrier()
		
	#########################################
	#	Print out time and finished jobs	#
	#########################################
	if rank == 0:
		remaining_time = str(datetime.timedelta(seconds=0)).split(".")[0]
		print(f"\rFinished Jobs: [{max_jobs}/{max_jobs}], Remaining time {remaining_time}            ", end='', flush=False)
		total_time = time.time() - start_time	# Total time
		total_time = str(datetime.timedelta(seconds=total_time)).split(".")[0] # Remove ms
		print(f"\n[NOTE] Finished {max_jobs} Jobs in {total_time}.")

	########################################################
	#	Read all the results and put it into npy files	#
	########################################################
	os.chdir(path)
	if rank == 0:
		start = time.time()
		print("[STATUS] Gathering results...", end='', flush=False)

		# Redefine tasks as now all the tasks are read
		tasks = create_task_folder_list(Map) # Structure tasks

		# Create shapes of the arrays which are filled and saved later
		log_tau, _,_,_,_,_,_,_,_,_,_ = sir.read_model(f"{tasks['folders'][0]}/best1.mod")
		#stokes_inv	= np.ones (shape=(Map[1]-Map[0]+1,Map[3]-Map[2]+1,4 ,data.shape[3]))
		stokes_inv = p.Profile(stk.nx, stk.ny, stk.nw)
		stokes_inv.wave = stk.wave # Copy wavelength positions
		models_inv1		= m.Model(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, npar=len(log_tau))
		models_inv2		= m.Model(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, npar=len(log_tau))
		errors_inv1		= m.Error(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, npar=len(log_tau))
		errors_inv2		= m.Error(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, npar=len(log_tau))
		best_guesses1	= m.Model(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, npar=len(log_tau))
		best_guesses2	= m.Model(nx = Map[1]-Map[0]+1, ny = Map[3]-Map[2]+1, npar=len(log_tau))
		chi2			= np.zeros(shape=(Map[1]-Map[0]+1,Map[3]-Map[2]+1))


		models_inv1.read_results(tasks, 'best1.mod', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		models_inv2.read_results(tasks, 'best2.mod', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		errors_inv1.read_results(tasks, 'best1.err', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		errors_inv2.read_results(tasks, 'best2.err', path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		best_guesses1.read_results(tasks, d.best_guess1, path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		best_guesses2.read_results(tasks, d.best_guess2, path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)
		stokes_inv.read_results(tasks, f"best.per", path, Map[1]-Map[0]+1, Map[3]-Map[2]+1)

		# Collect data from task folders and delete the folder
		for i in range(len(tasks['x'])):
			# Get folder name and x,y positions
			folder = tasks['folders'][i]
			x = tasks['x'][i]
			y = tasks['y'][i]
			
			# Read chi2 file
			chi2[x-Map[0],y-Map[2]] = sir.read_chi2(f"{folder}/{d.inv_trol_file[:d.inv_trol_file.rfind('.')]}.chi")

			# Remove folder
			shutil.rmtree(folder)

		# Create directory if / exists in inv_out in config
		if '/' in conf['inv_out']:
			temp = os.path.join(path, conf["inv_out"])
			if not exists(temp[:temp.rfind('/')]):
				os.mkdir(temp[:temp.rfind('/')])

		# Save as npy files
		#np.save(conf['inv_out'] + d.end_stokes, stokes_inv)
		stokes_inv.write(os.path.join(path,conf['inv_out']) + d.end_stokes)
		models_inv1.write(conf['inv_out'] + d.end_models1)
		models_inv2.write(conf['inv_out'] + d.end_models2)
		errors_inv1.write(conf['inv_out'] + d.end_errors1)
		errors_inv2.write(conf['inv_out'] + d.end_errors2)
		best_guesses1.write(d.best_guess1.replace(".mod",".bin"))
		best_guesses2.write(d.best_guess2.replace(".mod",".bin"))
		np.save(conf['chi2'],chi2)
		
		# Print needed time
		end = time.time()
		Time = str(datetime.timedelta(seconds=end-start)).split(".")[0]
		print(f" Finished in {Time}")

	comm.barrier()

	return chi2s, chi2s_num

