import sys
import os
import shutil
from os.path import exists
import sir
import definitions as d

import preprocess.merge
import preprocess.correction_spectral_veil

import _1C.inversion

import _2C.inversion

import _MC.create_models  # Creating Models
import _MC.add_noise	  # Adding Noise
import _MC.synthesis	  # Synthesis
import _MC.inversion	  # Inversion


from mpi4py import MPI

def initial(mode):
	"""
	Initial print outs and preparation

	Parameter
	---------
	mode : string
		Mode which is used
	
	Return
	------
	None
	"""
	print()
	print("╭───────────────────────────────────────────────────╮")
	print("│ SIR - MIG                                         │")
	print("│ Version 1.0                                       │")
	print("│ Multiple Initial Guesses                          │")
	print(f"│ Mode: {mode}                                          │")
	print("│ Author: Mike Moser                                │")
	print("╰───────────────────────────────────────────────────╯")	
	print()

if "-h" in sys.argv:
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
	sir.option("--no-inv","Do not perform the inversion [Only in mode MC]")
	sir.option("--only-inv","Only perform inversion [Only in mode MC]")
	sys.exit()


# Read the config file from the input
conf = sir.read_config(sys.argv[1])

# Implement MPI stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
	initial(conf['mode'])


#####################
# MERGING THE DATA	#
#####################
if (conf['mode'] == "1C" or conf['mode'] == "2C") and rank == 0:
	
	if conf['preprocess'] == "1":
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
				preprocess.merge.merge(conf, Dir, ends)
		######################
		# NORMALISE THE DATA #
		######################
		preprocess.normalise.normalise(conf)
		#########################
		# CORRECT SPECTRAL VEIL #
		#########################
		if conf['fts_file'] != '':
			preprocess.correction_spectral_veil.correct_spectral_veil(conf)
		elif not exists(os.path.join(conf['path'],conf['cube_inv'])) and len(conf['quiet_sun']) > 1:
			shutil.copy(os.path.join(conf['path'],conf['cube']) + d.end_norm, os.path.join(conf['path'],conf['cube_inv']))
		else:
			print("[main] No preprocessing steps are actually performed. Is 'preprocess=1' really needed in the config file?")

comm.barrier()
	
if conf['mode'] == "1C":
	#####################
	# PERFORM INVERSION #
	#####################
	_1C.inversion.inversion(conf, comm, rank, size)

if conf['mode'] == '2C':
	#####################
	# PERFORM INVERSION #
	#####################
	_2C.inversion.inversion(conf, comm, rank, size)

if conf['mode'] == 'MC':
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
				_MC.create_models.create_models(conf)
			print("[STATUS] Perform synthesis")

		comm.barrier()

		#####################
		# Perform Synthesis #
		#####################
		if not "--no-syn" in sys.argv:
			_MC.synthesis.synthesis(conf, comm, rank, size)

		comm.barrier()

		#############################
		# Add Noise to the Profiles #
		#############################
		if rank == 0:
			if "--no-noise" in sys.argv:
				print("-------> Skip adding noise ...")
			else:
				print("[STATUS] Add noise")
				_MC.add_noise.add_noise(conf, False)
			if not "--no-inv" in sys.argv:
				print("[STATUS] Perform inversion")
	else:
		if rank == 0:
			print("[STATUS] Only perform inversion")

	if not "--no-inv" in sys.argv:
		#####################
		# PERFORM INVERSION #
		#####################
		_MC.inversion.inversion(conf, comm, rank, size)




comm.barrier()
	
MPI.Finalize()

