#!/usr/bin/env python
"""
Perform synthesis of many profiles with simple parallelization using mpi4py
"""

import os, shutil
import glob
import sys
import sir
import model as m
import profile_stk as p
import definitions as d
import numpy as np

def help():
	"""
	Help page
	"""	
	print("synthesis - Executes the synthesis")
	print("Usage: python synthesis.py [OPTION]")
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


def synthesis(conf, comm, rank, size):
	"""
	Performs the synthesis of all the models.

	Parameter
	---------
	config : dict
		Dcitionary with config information
	rank, comm, size : MPI variables

	Return
	------
	None
	"""
	comm.barrier()

	###################################################
	#		READ PARAMETERS FROM CONFIG FILE		#	
	###################################################
	# Define parameters for easier understanding
	path = conf["path"]
	model = conf["model"]
	model1 = model.replace('.mod','') # For simplicity reasons defined

	cycles = conf['cycles'] # How many cycles
	line_file = conf['line'] # line file
	abundance_file = conf['abundance'] # Abundance file	
	end='.per' # ending of the file for the profile
	
	models = m.read_model(os.path.join(path, conf['model_out']))

	####################################
	#	CREATE GRID AND CONFIG FILE	#
	####################################
	if rank == 0:
		sir.write_control_mc(os.path.join(path,d.syn_trol_file), conf, Type="syn")
		sir.write_grid_mc(conf,os.path.join(path,d.Grid))
			
	if rank == 0:
		# Check if there are old task folders and delete them => can result to errors
		if len(glob.glob(os.path.join(path,d.task_start + "*"))) > 0:
			print("-------> Deleting old task folders")
			for i in glob.glob(os.path.join(path,d.task_start + "*")):
				shutil.rmtree(i)
	
	comm.barrier() # Continues when all processes entered this function

	
	####################################
	#		START INVERSION PART	#
	####################################
	num = conf['num']

	performed_models = 0 # Counts how many models are performed
	finished_jobs = 0

	tasks = create_task_folder_list(conf['num'])

	for i in range(rank, conf['num'], size):
		# Create task folder
		task_folder = os.path.join(path, tasks['folders'][i])
		os.makedirs(task_folder, exist_ok=True)
		
		# Copy SIR files and guess model
		sir_files = [d.syn_trol_file, "sir.x", conf['line'], d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(task_folder, sir_file))
		
		# Extract model from the npy file
		models.write_model(os.path.join(task_folder,d.model_syn), d.header, i, 0)
	
		# Perform synthesis
		os.chdir(task_folder)
		os.system("echo " + d.syn_trol_file + " | " + " ./sir.x >/dev/null 2>/dev/null")
		
		performed_models += 1

		if finished_jobs < (conf['num'] - conf['num'] % size):
			finished_jobs = comm.allreduce(performed_models, op=MPI.SUM)

		if rank == 0:
			print(f"\rTotal finished Jobs: {finished_jobs}", end='', flush=False)

		os.chdir('../') # Go back in case relative paths are used
	os.chdir(path)
	# Collect data and save it
	if rank == 0:
		print(f"\rTotal finished Jobs: {conf['num']}", end='', flush=False)
		output =  os.path.join(path,conf["syn_out"]) # Output file
		atoms = [i.split(",") for i in conf['atoms']]

		# Read the profiles
		stk = p.Profile(conf['num'],1,0)
		stk.read_results_MC(path, tasks, d.profile)
		
		stk.write(f"{conf['syn_out']}")
		
		for i in range(conf['num']):
			shutil.rmtree(tasks['folders'][i])
		print(f"\r-------> Finished with {conf['num']} synthesised models.")

	comm.barrier()
