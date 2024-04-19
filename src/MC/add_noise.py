"""
Adding noise to many synthesis results
"""

import numpy as np 
import sys
import os
import sir


def help() -> None:
	"""
	Help Page
	"""
	print("add_noise - Adding noise to synthesis results of a MC simulation")
	print("Usage: python add_noise.py [OPTION]")
	print()
	sir.option("[1. Pos]","Config file")
	sys.exit()


def add_noise(conf: dict, verbose: bool = True) -> None:
	"""
	Adds noise to the Stokes Profiles

	Parameter
	---------
	conf : dict
		Information from the config file
	verbose : bool, optional
		Print out status. Default: True

	Return
	------
	None


	"""
	if verbose:
		print("[STATUS] Add noise to the synthesised profiles")
	##################
	# READ ARGUMENTS #
	##################
	path = conf["path"]
	Input = os.path.join(path, conf["syn_out"])			# Input synthesised model profiles/syn
	Output = os.path.join(path, conf["noise_out"])		# Generic output profiles/noise

	noise_I = conf['noise_I']  # Noise in I
	noise_Q = conf['noise_Q']  # Noise in Q
	noise_U = conf['noise_U']  # Noise in U
	noise_V = conf['noise_V']  # Noise in V

	atoms = [i.split(",") for i in conf['atoms']]

	#############################
	# ADD NOISE TO EACH PROFILE #
	#############################
	syn = np.load(f"{Input}")

	noise_Is = np.random.normal(scale=float(noise_I), size=(syn.shape[0], syn.shape[2]))
	noise_Qs = np.random.normal(scale=float(noise_Q), size=(syn.shape[0], syn.shape[2]))
	noise_Us = np.random.normal(scale=float(noise_U), size=(syn.shape[0], syn.shape[2]))
	noise_Vs = np.random.normal(scale=float(noise_V), size=(syn.shape[0], syn.shape[2]))

	syn[:, 2, :] = syn[:, 2, :] + noise_Is
	syn[:, 3, :] = syn[:, 3, :] + noise_Qs
	syn[:, 4, :] = syn[:, 4, :] + noise_Us
	syn[:, 5, :] = syn[:, 5, :] + noise_Vs

	np.save(f"{Output}", syn)

	return


if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	if len(sys.argv) < 2:
		print("[ERROR] No config file provided!")
		sys.exit()
	elif not os.path.exists(sys.argv[1]):
		print("[ERROR] Config file does not exist!")
		sys.exit()
	add_noise(sir.read_config(sys.argv[1]))
