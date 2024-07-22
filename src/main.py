#!/usr/bin/env python3
import sys
import os
from os.path import exists
import sir

def help():
	"""
	Help Page

	Parameters
	----------
	None

	Returns
	-------
	None
	"""
	print("main - Executes all the scripts (merging, normalisation, spectral veil correction, inversion)")
	print("       Normalisation and spectral veil correciton is only performed if preprocess is active in the config file!")
	print("Usage: mpirun -np X python inversion.py [OPTION]")
	print()
	sir.option("1:","Config file")
	sir.option("-dir:","Merge data with the files in this directory (if not set, data is not merged) [Only in mode 1C and 2C]")
	sir.option("-end:","Ending for GRIS data (001,002,...) [Only in mode 1C and 2C]")
	sir.option("-save:","Additional savepath for plots from spectral veil correction, ending with / [Only in mode 1C and 2C]")
	sir.option("--no-create","Do not create models [Only in mode MC]")
	sir.option("--no-syn","Do not perform the synthesis [Only in mode MC]")
	sir.option("--no-noise","Do not add noise [Only in mode MC]")
	sir.option("--no-inv","Do not perform the inversion")
	sir.option("--only-inv","Only perform inversion [Only in mode MC]")
	sir.option("--debug", "Debugging (task folders in inversions are not deleted)")
	sir.option("--no-progress", "Do not print a progress bar")
	sys.exit()

def sir_mig():
	"""
	Function which executes the programme

	Parameters
	----------
	None

	Returns
	-------
	None
	"""
	# Implement MPI stuff
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()

	# Read the config file from the input
	if rank == 0:
		conf = sir.read_config(sys.argv[1])
	else:
		conf = None
	conf = comm.bcast(conf, root=0)

	if rank == 0:
		sir.initial(conf['mode'])
		if("--debug" in sys.argv):
			print("-------> Flag 'debug' used")
			
	comm.barrier()

	#####################
	# PREPROCESS	#
	#####################
	if (conf['mode'] == "1C" or conf['mode'] == "2C") and rank == 0:
		# Check for flags and print it out
		if "--no-inv" in sys.argv:
			print("-------> Flag 'no-inv' used")
		if conf['preprocess'] == "1":
			import preprocess
			print("[STATUS] Preprocess data")
			if exists(os.path.join(conf['path'],conf['cube'])):
					print("-------> Skipping preprocessing data (already performed)")
			else:
				if "-dir" in sys.argv and conf['preprocess'] == "1":
					print("[STATUS] Merge the data")
					Dir = sys.argv[sys.argv.index("-dir")+1]
					
					pro = preprocess.merge(Dir, conf["ending"], conf["instrument"], conf["path"], conf["shift_wave"], conf["save_cube"] == "1")

				######################
				# NORMALISE THE DATA #
				######################
				pro = preprocess.normalise(pro, conf["instrument"], conf["quiet_sun"], conf["path"], conf["save_cube"] == "1")
				
				#########################
				# CORRECT SPECTRAL VEIL #
				#########################
				if conf['fts_file'] != '':
					preprocess.correct_spectral_veil_conf(pro, conf)
				elif not exists(os.path.join(conf['path'],conf['cube'])):
					print("-------> Saving data (this might take a while) ...")
					pro.write(os.path.join(conf['path'],conf['cube']))
				else:
					print("[main] No preprocessing steps are actually performed. Is 'preprocess=1' really needed in the config file?")

	comm.barrier()
		
	if conf['mode'] == "1C":
		#####################
		# PERFORM INVERSION #
		#####################
		if not "--no-inv" in sys.argv:
			import inversion
			inversion.inversion_1c(conf, comm, rank, size, MPI, "--debug" in sys.argv, not "--no-progress" in sys.argv)

	elif conf['mode'] == '2C':
		#####################
		# PERFORM INVERSION #
		#####################
		if not "--no-inv" in sys.argv:
			import inversion
			inversion.inversion_2c(conf, comm, rank, size, MPI, "--debug" in sys.argv, not "--no-progress" in sys.argv)

	elif conf['mode'] == 'MC':
		import simulation  # Creating Models
		import inversion
		if rank == 0:
			# Check for flags and print it out
			if "--no-create" in sys.argv:
				print("-------> Flag 'no-create' used")
			if "--no-syn" in sys.argv:
				print("-------> Flag 'no-syn' used")
			if "--no-noise" in sys.argv:
				print("-------> Flag 'no-noise' used")
			if "--no-inv" in sys.argv:
				print("-------> Flag 'no-inv' used")
			if "--only-inv" in sys.argv:
				print("-------> Flag 'only-inv' used")
			print("[STATUS] Start Simulation")
		if not "--only-inv" in sys.argv:
			#################
			# Create Models #
			#################
			if rank == 0:
				if "--no-create" in sys.argv:
					print("-------> Skip creating models ...")
				else:
					print("[STATUS] Create Models")
					simulation.create_models(conf)
				print("[STATUS] Perform synthesis")

			comm.barrier()

			#####################
			# Perform Synthesis #
			#####################
			if not "--no-syn" in sys.argv:
				simulation.synthesis(conf, comm, rank, size, MPI, not "--no-progress" in sys.argv)

			comm.barrier()

			#############################
			# Add Noise to the Profiles #
			#############################
			if rank == 0:
				if "--no-noise" in sys.argv:
					print("-------> Skip adding noise ...")
				else:
					print("[STATUS] Add noise")
					simulation.add_noise(conf, False)
		else:
			if rank == 0:
				print("[STATUS] Only perform inversion")

		if not "--no-inv" in sys.argv:
			#####################
			# PERFORM INVERSION #
			#####################
			inversion.inversion_mc(conf, comm, rank, size, MPI, "--debug" in sys.argv, not "--no-progress" in sys.argv)




	comm.barrier()
		
	MPI.Finalize()

if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	else:
		sir_mig()


