"""
Cleans the files from the inversion
"""

import sys
import os
sys.path.append(sys.path[0] + "/../")
import sir
import definitions as d



if "-h" in sys.argv or (len(sys.argv) < 2):
	print("clean - Cleans up the inversion")
	print("Usage: python clean.py [OPTION]")
	print()
	sir.option("[1. Pos]","Config file")
	sys.exit()

#######################################################################
#	    READ INPUT, WRITE DEFINITIONS AND LOAD DATA			  #
#######################################################################
conf    = sir.read_config(sys.argv[1])
path    = conf["path"]
inv_out = os.path.join(path,conf["inv_out"])
if conf['mode'] == "1C":
	os.system(f"rm -rv {inv_out+d.end_stokes}")
	os.system(f"rm -rv {inv_out+d.end_models}")
	os.system(f"rm -rv {inv_out+d.end_errors}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess_file)}")
	os.system(f"rm -rv {os.path.join(path,conf['chi2'])}")
	if conf['psf'] != '':
		os.system(f"rm -rv {os.path.join(path,d.psf)}")

elif conf['mode'] == '2C':
	os.system(f"rm -rv {inv_out+d.end_stokes}")
	os.system(f"rm -rv {inv_out+d.end_models1}")
	os.system(f"rm -rv {inv_out+d.end_models2}")
	os.system(f"rm -rv {inv_out+d.end_errors1}")
	os.system(f"rm -rv {inv_out+d.end_errors2}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess1_file)}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess2_file)}")
	os.system(f"rm -rv {os.path.join(path,conf['chi2'])}")
	if conf['psf'] != '':
		os.system(f"rm -rv {os.path.join(path,d.psf)}")

elif conf['mode'] == 'MC':
	inv_out = os.path.join(path,conf["inv_out"])
	guess_out = os.path.join(path,d.best_guess_file)
	os.system(f"rm -rv {inv_out+d.end_stokes}")
	os.system(f"rm -rv {inv_out+d.end_models}")
	os.system(f"rm -rv {inv_out+d.end_errors}")
	os.system(f"rm -rv {os.path.join(path,d.best_guess_file)}")
	os.system(f"rm -rv {os.path.join(path,conf['chi2'])}")
	os.system(f"rm -rf {guess_out}")

else:
	print("[clean] Mode not known")
