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
	sir.option("-ending:","Ending for GRIS data (001,002,...) [Only in mode 1C and 2C]")
	sir.option("-save:","Additional savepath for plots from spectral veil correction, ending with / [Only in mode 1C and 2C]")
	sir.option("--no-create","Do not create models [Only in mode MC]")
	sir.option("--no-syn","Do not perform the synthesis [Only in mode MC]")
	sir.option("--no-noise","Do not add noise [Only in mode MC]")
	sir.option("--no-inv","Do not perform the inversion")
	sir.option("--only-inv","Only perform inversion [Only in mode MC]")
	sys.exit()

def main():
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
	conf = sir.read_config(sys.argv[1])
	if rank == 0:
		sir.initial(conf['mode'])

	comm.barrier()

	#####################
	# PREPROCESS	#
	#####################
	if (conf['mode'] == "1C" or conf['mode'] == "2C") and rank == 0:
		
		if conf['preprocess'] == "1":
			import preprocess
			print("[STATUS] Preprocess data")
			if exists(os.path.join(conf['path'],conf['cube_inv'])):
					print("-------> Skipping preprocessing data (already performed)")
			else:
				if "-dir" in sys.argv and conf['preprocess'] == "1":
					print("[STATUS] Merge the data")
					Dir = sys.argv[sys.argv.index("-dir")+1]
					if "-ending" in sys.argv:
						ends = sys.argv[sys.argv.index("-ending")+1]
					elif conf['instrument'] == 'GRIS':
						ends = "001"
					else:
						ends = ''
					preprocess.merge(conf, Dir, ends)
			######################
			# NORMALISE THE DATA #
			######################
			preprocess.normalise(conf)
			#########################
			# CORRECT SPECTRAL VEIL #
			#########################
			if conf['fts_file'] != '':
				preprocess.correct_spectral_veil(conf)
			elif not exists(os.path.join(conf['path'],conf['cube_inv'])) and len(conf['quiet_sun']) > 1:
				import shutil
				import definitions as d
				shutil.copy(os.path.join(conf['path'],conf['cube']).replace(".bin","") + d.end_norm, os.path.join(conf['path'],conf['cube_inv']))
			else:
				print("[main] No preprocessing steps are actually performed. Is 'preprocess=1' really needed in the config file?")

	comm.barrier()
		
	if conf['mode'] == "1C":
		#####################
		# PERFORM INVERSION #
		#####################
		if not "--no-inv" in sys.argv:
			import inversion
			inversion.inversion_1c(conf, comm, rank, size, MPI)

	elif conf['mode'] == '2C':
		#####################
		# PERFORM INVERSION #
		#####################
		if not "--no-inv" in sys.argv:
			import inversion
			inversion.inversion_2c(conf, comm, rank, size, MPI)

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
				simulation.synthesis(conf, comm, rank, size, MPI)

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
			inversion.inversion_mc(conf, comm, rank, size, MPI)




	comm.barrier()
		
	MPI.Finalize()

if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	else:
		main()

