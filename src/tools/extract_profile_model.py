"""
Extracts the profile file and the mdoel file from a data cube
"""
import numpy as np
import sys, os, shutil
sys.path.append(sys.path[0]+"/..")
import sir
import definitions as d
from os.path import exists
import profile_stk as p
import model_atm as m

def help():
	"""
	Prints information how to use this script.
	"""
	print("extract_profile_model - Extracts the profile file and the model file from a data cube")
	print("Usage: python extract_profile_model [OPTION]")
	print()
	sir.option("[1. Pos.]","Config")
	sir.option("[2. Pos]","x position in absolute position")
	sir.option("[3. Pos.]","y position in absolute position [only for mode 1C and 2C]")
	print()
	sir.option("-save","Additional save path (optional, default './')")
	sir.option("-add","Additional text in filenames (optional)")

	sys.exit()


def extract_profile_model_1C(conf, x, y):
	"""
	Extract the profiles and the models from a specific pixel position.
	It saves and copies all the necessary files to rerun the inversion.
	
	Parameter
	---------
	config : dict
		Dictionary containing all the information from the config
	x : int
		x-position of the pixel in absolute position (taking into account the full observed map) 
	y : int
		y-position of the pixel in absolute position (taking into account the full observed map)

	Return
	------
	None

	"""

	path = conf["path"]

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath):
			os.mkdir(savepath)

	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Load wavelenghts and map of (x,y) used in the inversion
	Map = conf['map']
	if x < Map[0] or x > Map[1]:
		print("[ERROR] x is not in the range (%i,%i)!" % (Map[0],Map[1]))
		sys.exit(1)
	if y < Map[2] or y > Map[3]:
		print("[ERROR] y is not in the range (%i,%i)!" % (Map[2],Map[3]))
		sys.exit(1)

	# Load data
	obs1 = p.Profile(os.path.join(path, conf['cube_inv']))
	obs1.cut_to_wave(conf["range_wave"]) # Cut the values to data used in the inversion

	inv = p.Profile(os.path.join(path, conf['inv_out']) + d.end_stokes)
	mod = m.Model(filename=os.path.join(path,conf['inv_out'] + d.end_models))
	guess = m.Model(filename=os.path.join(path,d.best_guess.replace(".mod",".bin")))
	err = m.Model(filename=os.path.join(path,conf['inv_out'] + d.end_errors))

	# Change to the reduced Map
	x = x - Map[0]
	y = y - Map[2]
	
	# Save data in formats for SIR
	obs1.write_profile(savepath + "profile" + add + ".per", x + Map[0], y + Map[2], d.Grid)
	inv.write_profile(savepath + "profile_result" + add + ".per", x, y, d.Grid)
	mod.write_model(savepath + "model_result" + add + ".mod", x, y)
	guess.write_model(savepath + conf["model"].replace(".mod","") + add + ".mod", x, y)
	err.write_model(savepath + "model_result" + add + ".err", x, y)
	
	# Copy stuff for the inversion
	if savepath != '':
		if conf['psf'] != '':
			sir_files = [d.inv_trol_file, "sir.x", conf['line'], d.Grid, conf['abundance'], conf['psf']]
		else:
			sir_files = [d.inv_trol_file, "sir.x", conf['line'], d.Grid, conf['abundance']]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(savepath, sir_file))

def extract_profile_model_MC(conf, num):
	"""
	Extract the profiles and the models from a specific pixel position.
	It saves and copies all the necessary files to rerun the inversion.
	
	Parameters
	----------
	config : dict
		Dictionary containing all the information from the config
	num : int
		Model number

	Returns
	-------
	None

	"""

	path = conf["path"]
	model = conf['model']
	line_file = conf['line']	
	abundance_file = conf['abundance']

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath):
			os.mkdir(savepath)

	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Load data
	noise = p.Profile(os.path.join(path,conf['noise_out']))
	syn = p.Profile(os.path.join(path,conf['syn_out']))
	inv = p.Profile(os.path.join(path,f"{conf['inv_out']}{d.end_stokes}"))

	noise.write_profile_mc(os.path.join(savepath,'noise.per' + add), num-1)
	syn.write_profile_mc(os.path.join(savepath,'syn.per' + add), num-1)
	inv.write_profile_mc(os.path.join(savepath,'inv.per' + add), num-1)
	

	# Now with the class model/error
	mod = m.Model(os.path.join(path,conf['inv_out'] + d.inv_models))
	gue = m.Model(os.path.join(path,d.best_guess.replace('.mod','.bin')))
	syn = m.Model(os.path.join(path,conf['model_out']))
	err = m.Model(os.path.join(path, os.path.join(path,conf['inv_out'] + d.inv_errors)))
	
	syn.write_model(os.path.join(savepath,'syn.mod' + add), num-1,0)
	mod.write_model(os.path.join(savepath,'res.mod' + add), num-1,0)
	gue.write_model(os.path.join(savepath,'guess.mod' + add),num-1,0)
	err.write_model(os.path.join(savepath,'res.err' + add), num-1,0)

	#sir.write_model_npy(os.path.join(savepath,'syn.mod' + add), d.header,syn[num-1,:,:]) # Old


	# Copy stuff for the inversion
	if savepath != '':
		sir_files = [d.inv_trol_file, model, "sir.x", line_file, d.Grid, abundance_file]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(savepath, sir_file))

def extract_profile_model_2C(conf, x, y):
	"""
	Extract the profiles and the models from a specific pixel position.
	It saves and copies all the necessary files to rerun the inversion.
	
	Parameter
	---------
	conf : dict
		Dictionary containing all the information from the config
	x : int
		x-position of the pixel in absolute position (taking into account the full observed map) 
	y : int
		y-position of the pixel in absolute position (taking into account the full observed map)

	Return
	------
	None

	"""
	path = conf["path"]

	# Additional savepath
	savepath = ''
	if '-save' in sys.argv:
		savepath = sys.argv[sys.argv.index("-save")+1]
		if not exists(savepath):
			os.mkdir(savepath)

	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	Map = conf['map']
	
	# Load data
	obs1 = p.read_profile(conf, filename=conf['cube_inv'])	
	inv = p.read_profile(os.path.join(path,conf['inv_out'] + d.end_stokes))
	guess1 = m.read_model(os.path.join(path,d.best_guess1.replace(".mod",".bin")))
	guess2 = m.read_model(os.path.join(path,d.best_guess2.replace(".mod",".bin")))
	mod1 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models1))
	mod2 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models2))
	err1 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_errors1))
	err2 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_errors2))
	

	# Cut observation
	obs1.cut_to_map(Map)
	obs1.cut_to_wave(conf["range_wave"]) # Cut the values to data used in the inversion

	# Change to the reduced Map
	x = x - Map[0]
	y = y - Map[2]
	
	# Save data in formats for SIR
	obs1.write_profile(savepath + "profile" + add + ".per", x, y,d.Grid)
	inv.write_profile(savepath + "profile_result" + add + ".per", x, y)
	guess1.write_model(savepath + "guess1" + add + ".mod",x,y)
	guess2.write_model(savepath + "guess2" + add + ".mod",x,y)
	mod1.write_model(savepath + "model_result1" + add + ".mod", x,y)
	mod2.write_model(savepath + "model_result2" + add + ".mod", x,y)
	err1.write_model(savepath + "model_result1" + add + ".err", x,y)
	err2.write_model(savepath + "model_result2" + add + ".err", x,y)

	# Copy stuff for the inversion
	if savepath != '':
		if conf['psf'] != '':
			sir_files = [d.inv_trol_file, "sir.x", conf['line'], d.Grid, conf['abundance'], conf['psf']]
		else:
			sir_files = [d.inv_trol_file, "sir.x", conf['line'], d.Grid, conf['abundance']]
		for sir_file in sir_files:
			shutil.copy(os.path.join(path, sir_file), os.path.join(savepath, sir_file))

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1])

	if conf['mode'] == "MC":
		extract_profile_model_MC(conf, int(sys.argv[2]))
	elif conf['mode'] == "1C":
		extract_profile_model_1C(conf, int(sys.argv[2]), int(sys.argv[3]))
	elif conf['mode'] == "2C":
		extract_profile_model_2C(conf, int(sys.argv[2]), int(sys.argv[3]))



