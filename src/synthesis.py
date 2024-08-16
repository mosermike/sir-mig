
import definitions as d
import glob
import model_atm as m
import os
import profile_stk as p
import sir
import shutil

def synthesis(conf : dict, comm, rank : int, size : int, MPI, debug : bool=False, progress : bool = True):
	"""
	Performs the synthesis of all the models.

	Parameters
	----------
	config : dict
		Dcitionary with config information
	rank : int
		Number of this process
	comm : Intracomm
		Communicator from MPI
	size : int
		Number of processes
	MPI : library
		Library MPI
	debug : bool,optional
		Do not delete created files
	progress : bool,optional
		Print a progress bar

	Returns
	-------
	None

	"""
	
	comm.barrier()

	###################################################
	#		READ PARAMETERS FROM CONFIG FILE		#	
	###################################################
	# Define parameters for easier understanding
	path = conf["path"]
	abundance_file = conf['abundance'] # Abundance file	
	
	models = m.read_model(os.path.join(path, conf['syn_in']))

	####################################
	#	CREATE GRID AND CONFIG FILE	#
	####################################
	if rank == 0:
		print(f"-------> Write control file")
		sir.write_control(os.path.join(path,d.syn_trol_file), conf)
		print(f"-------> Write grid file")
		sir.write_grid(conf,os.path.join(path,d.Grid))
			
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

	performed_models = 0 # Counts how many models are performed
	finished_jobs = 0

	tasks = sir.create_task_folder_list(int(models.nx))

	for i in range(rank, models.nx, size):
		# Create task folder
		task_folder = os.path.join(path, tasks['folders'][i])
		os.makedirs(task_folder, exist_ok=True)
		
		# Copy SIR files and guess model
		sir_files = [d.syn_trol_file, "sir.x", conf['line'], d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(task_folder, sir_file))
		
		# Extract model
		models.write_model(os.path.join(task_folder,d.model_syn), i, 0)
	
		# Perform synthesis
		os.chdir(task_folder)
		os.system("echo " + d.syn_trol_file + " | " + " ./sir.x 2>1 >syn.log")
		
		performed_models += 1

		if finished_jobs < (models.nx - models.nx % size):
			finished_jobs = comm.allreduce(performed_models, op=MPI.SUM)

		if rank == 0 and progress:
			print(f"\rTotal finished Jobs: {finished_jobs}", end='', flush=False)

		os.chdir('../') # Go back in case relative paths are used
	
	# Collect data and save it
	if rank == 0:
		print(f"\rTotal finished Jobs: {models.nx}", end='', flush=False)

		# Read the profiles
		stk = p.profile_stk(models.nx,1,0)
		stk.read_results_MC(path, tasks, d.profile)
		
		stk.write(f"{os.path.join(conf['path'],conf['syn_out'])}")
		
		if not debug:
			for i in range(models.nx):
				shutil.rmtree(tasks['folders'][i])
			os.remove(os.path.join(path,d.syn_trol_file))
			os.remove(os.path.join(path,d.Grid))
		else:
			print("-------> No created files are deleted. Debug option active.")

		print(f"\r-------> Finished with {models.nx} synthesised models.")

	comm.barrier()