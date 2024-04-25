"""
Removes all created files
"""

import sys
import os
sys.path.append(sys.path[0]+"/..")
import sir
import definitions as d
# Import library


if "-h" in sys.argv or (len(sys.argv) < 2):
	print("distclean - Cleans up all created files")
	print("Usage: python distclean.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Config file")
	
	sys.exit()

#######################################################################
#	    READ INPUT, WRITE DEFINITIONS AND LOAD DATA			  #
#######################################################################
conf = sir.read_config(sys.argv[1])
path = conf["path"]
inv_out = os.path.join(path,conf["inv_out"])

# Check if path exists
if conf['path'] != os.path.abspath(os.getcwd()):
	Inp = input(f"[NOTE] Paths are not the same. You are in {os.path.abspath(os.getcwd())} but the config path is {conf['path']}. Do you want to continue? [y/n] ")
	if Inp != "y":
		print('Abort')
		sys.exit()

if conf['mode'] == "1C":

	os.system(f"rm -rv {os.path.join(path,conf['cube'].replace('.npy','')+d.end_norm + '.npy')}")
	os.system(f"rm -rv {inv_out+d.end_stokes}")
	os.system(f"rm -rv {inv_out+d.end_models}")
	os.system(f"rm -rv {inv_out+d.end_errors}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess.replace('.mod','.npy'))}")
	os.system(f"rm -rv {os.path.join(path,conf['chi2'])}")
	os.system(f"rm -rv {os.path.join(path,conf['waves'])}")
	os.system(f"rm -rv {os.path.join(path,d.inv_trol_file)}")
	os.system(f"rm -rv {os.path.join(path,d.Grid)}")
	os.system(f"rm -rv {os.path.join(path,d.veil_parameters)}")
	os.system(f"rm -rv {os.path.join(path,d.header_infos)}")

	if conf['psf'] != '':
		os.system(f"rm -rv {os.path.join(path,conf['psf'])}")

elif conf['mode'] == '2C':
	os.system(f"rm -rv {os.path.join(path,conf['cube'].replace('.npy','')+d.end_norm + '.npy')}")
	os.system(f"rm -rv {inv_out+d.end_stokes}")
	os.system(f"rm -rv {inv_out+d.end_models1}")
	os.system(f"rm -rv {inv_out+d.end_models2}")
	os.system(f"rm -rv {inv_out+d.end_errors1}")
	os.system(f"rm -rv {inv_out+d.end_errors2}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess1.replace('.mod','.npy'))}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess2.replace('.mod','.npy'))}")
	os.system(f"rm -rv {os.path.join(path,conf['chi2'])}")
	os.system(f"rm -rv {os.path.join(path,conf['waves'])}")
	os.system(f"rm -rv {os.path.join(path,d.inv_trol_file)}")
	os.system(f"rm -rv {os.path.join(path,d.Grid)}")
	os.system(f"rm -rv {os.path.join(path,d.veil_parameters)}")
	os.system(f"rm -rv {os.path.join(path,d.header_infos)}")

	if conf['psf'] != '':
		os.system(f"rm -rv {os.path.join(path,conf['psf'])}")

elif conf['mode'] == "MC":
	if conf["guess"] != '':
		guess   = os.path.join(path,conf["guess"])
		os.system(f"rm -rf {guess}")
	os.system(f"rm -rv {os.path.join(path,d.syn_trol_file)}")
	os.system(f"rm -rv {os.path.join(path,d.inv_trol_file)}")
	os.system(f"rm -rf {os.path.join(path, d.Grid)}")
	os.system(f"rm -rf {os.path.join(path, conf['chi2'])}")
	os.system(f"rm -rf {os.path.join(path, conf['model_out'])}")
	os.system(f"rm -rf {os.path.join(path,conf["syn_out"])}")
	os.system(f"rm -rf {os.path.join(path,conf["noise_out"])}")
	os.system(f"rm -rf {inv_out + d.end_stokes}*")
	os.system(f"rm -rf {inv_out + d.end_models}")
	os.system(f"rm -rf {inv_out + d.end_errors}")
	os.system(f"rm -rf {os.path.join(path,d.best_guess.replace('.mod','.npy'))}")

else:
	print("[distclean] Mode unknown")