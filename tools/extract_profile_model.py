"""
Extracts the profile file and the mdoel file from a data cube
"""
import numpy as np
import sys, os, shutil
sys.path.append(sys.path[0]+"/../src")
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


def extract_profile_model_1C(conf : dict, x : int, y : int, savepath : str):
	"""
	Extract the profiles and the models from a specific pixel position.
	It saves and copies all the necessary files to rerun the inversion.
	
	Parameters
	----------
	config : dict
		Dictionary containing all the information from the config
	x : int
		x-position of the pixel in absolute position (taking into account the full observed map) 
	y : int
		y-position of the pixel in absolute position (taking into account the full observed map)
	savepath : str
		Where the files are stored

			Return
	------
	None

	"""

	path = conf["path"]

	# Additional savepath
	if not exists(savepath[:savepath.rfind("/")]):
		os.mkdir(savepath[:savepath.rfind("/")])

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
	obs1 = p.read_profile(os.path.join(path, conf['cube']))
	obs1.cut_to_wave(conf["range_wave"]) # Cut the values to data used in the inversion
	obs1.cut_to_map(conf["map"]) # Cut the values to data used in the inversion

	inv   = p.read_profile(os.path.join(path, conf['inv_out']) + d.end_stokes)
	mod   = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models))
	guess = m.read_model(os.path.join(path,conf['inv_out'] + d.best_guess_file))
	err   = m.read_model(os.path.join(path,conf['inv_out'] + d.end_errors))

	inv.data_cut_wave = True
	inv._data_cut_map = True

	# Change to the reduced Map
	x = x - Map[0]
	y = y - Map[2]
	
	sir.write_grid(conf, os.path.join(savepath, d.Grid), obs1.wave)
	# Save data in formats for SIR
	obs1.write_profile(savepath + "profile" + add + ".per", x, y, os.path.join(savepath, d.Grid))
	inv.write_profile(savepath + "profile_result" + add + ".per", x, y, os.path.join(savepath, d.Grid))
	mod.write_model(savepath + "model_result" + add + ".mod", x, y)
	guess.write_model(savepath + d.guess.replace(".mod","") + add + ".mod", x, y)
	err.write_model(savepath + "model_result" + add + ".err", x, y)
	
	# Copy stuff for the inversion
	if conf['psf'] != '':
		if "gauss" in conf['psf']:
			sir.write_gauss_psf(float(conf['psf'].split(" ")[1]), os.path.join(savepath, d.psf))
		else:
			shutil.copy(os.path.join(path, conf['psf']), os.path.join(savepath, d.psf))
	
	sir_files = ["sir.x", conf['line'], conf['abundance']]
	for sir_file in sir_files:
		shutil.copy(os.path.join(path, sir_file), os.path.join(savepath, sir_file))
	
	sir.write_control(os.path.join(savepath, d.inv_trol_file), conf)

def extract_profile_model_MC(conf : dict, num : int, savepath : str):
	"""
	Extract the profiles and the models from a specific pixel position.
	It saves and copies all the necessary files to rerun the inversion.
	
	Parameters
	----------
	config : dict
		Dictionary containing all the information from the config
	num : int
		Model number
	savepath : str
		Where the files are stored

	Returns
	-------
	None

	"""

	path = conf["path"]
	model = conf['model']
	line_file = conf['line']	
	abundance_file = conf['abundance']

	# Additional savepath
	if not exists(savepath[:savepath.rfind("/")]):
		os.mkdir(savepath[:savepath.rfind("/")])

	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	# Load data
	noise = p.read_profile(os.path.join(path,conf['noise_out']+ d.end_stokes))
	syn = p.read_profile(os.path.join(path,conf["syn_out"]+ d.end_stokes))
	inv = p.read_profile(os.path.join(path,f"{conf['inv_out']}{d.end_stokes}"))

	noise.data_cut_wave = True
	noise._data_cut_map = True
	syn.data_cut_wave = True
	syn._data_cut_map = True
	inv.data_cut_wave = True
	inv._data_cut_map = True

	noise.write_profile_mc(os.path.join(savepath,'noise.per' + add), num-1)
	syn.write_profile_mc(os.path.join(savepath,'syn.per' + add), num-1)
	inv.write_profile_mc(os.path.join(savepath,'inv.per' + add), num-1)
	

	# Now with the class model/error
	mod = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models))
	gue = m.read_model(os.path.join(path,conf['inv_out'] + d.best_guess_file))
	syn = m.read_model(os.path.join(path,conf["syn_out"]+ d.end_models))
	err = m.read_model(os.path.join(path, os.path.join(path,conf['inv_out'] + d.end_errors)))
	
	syn.write_model(os.path.join(savepath,'syn.mod' + add), num-1,0)
	mod.write_model(os.path.join(savepath,'res.mod' + add), num-1,0)
	gue.write_model(os.path.join(savepath,'guess.mod' + add),num-1,0)
	err.write_model(os.path.join(savepath,'res.err' + add), num-1,0)

	#sir.write_model_npy(os.path.join(savepath,'syn.mod' + add), d.header,syn[num-1,:,:]) # Old


	# Copy stuff for the inversion
	sir_files = [model, "sir.x", line_file, abundance_file]
	for sir_file in sir_files:
		shutil.copy(os.path.join(path, sir_file), os.path.join(savepath, sir_file))

	sir.write_control(os.path.join(savepath, d.syn_trol_file), conf, "syn")
	sir.write_control(os.path.join(savepath, d.inv_trol_file), conf, "inv")
	sir.write_grid(conf, os.path.join(savepath, d.Grid))

def extract_profile_model_2C(conf : dict, x : int, y : int, savepath : str):
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
	savepath : str
		Where the files are stored

	Return
	------
	None

	"""
	path = conf["path"]

	# Additional savepath
	if not exists(savepath[:savepath.rfind("/")]):
		os.mkdir(savepath[:savepath.rfind("/")])

	# Additional text
	add = ''
	if '-add' in sys.argv:
		add = sys.argv[sys.argv.index("-add")+1]

	Map = conf['map']
	
	# Load data
	obs1 = p.read_profile(os.path.join(path,conf['cube']))
	inv = p.read_profile(os.path.join(path,conf['inv_out'] + d.end_stokes))
	guess1 = m.read_model(os.path.join(path,conf['inv_out'] + d.best_guess1_file))
	guess2 = m.read_model(os.path.join(path,conf['inv_out'] + d.best_guess2_file))
	mod1 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models1))
	mod2 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_models2))
	err1 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_errors1))
	err2 = m.read_model(os.path.join(path,conf['inv_out'] + d.end_errors2))
	

	# Cut observation
	obs1.cut_to_map(Map)
	obs1.cut_to_wave(conf["range_wave"]) # Cut the values to data used in the inversion
	inv.data_cut_wave = True
	inv._data_cut_map = True
	
	# Change to the reduced Map
	x = x - Map[0]
	y = y - Map[2]
	
	sir.write_grid(conf, os.path.join(savepath, d.Grid), obs1.wave)

	# Save data in formats for SIR
	obs1.write_profile(savepath + "profile" + add + ".per", x, y,os.path.join(savepath, d.Grid))
	inv.write_profile(savepath + "profile_result" + add + ".per", x, y, os.path.join(savepath, d.Grid))
	guess1.write_model(savepath + "guess1" + add + ".mod",x,y)
	guess2.write_model(savepath + "guess2" + add + ".mod",x,y)
	mod1.write_model(savepath + "model_result1" + add + ".mod", x,y)
	mod2.write_model(savepath + "model_result2" + add + ".mod", x,y)
	err1.write_model(savepath + "model_result1" + add + ".err", x,y)
	err2.write_model(savepath + "model_result2" + add + ".err", x,y)

	# Copy stuff for the inversion
	if conf['psf'] != '':
		if "gauss" in conf['psf']:
			sir.write_gauss_psf(float(conf['psf'].split(" ")[1]), os.path.join(savepath, d.psf))
		else:
			shutil.copy(os.path.join(path, conf['psf']), os.path.join(savepath, d.psf))

	sir_files = ["sir.x", conf['line'], conf['abundance']]
	
	for sir_file in sir_files:
		shutil.copy(os.path.join(path, sir_file), os.path.join(savepath, sir_file))
	
	sir.write_control(os.path.join(savepath, d.inv_trol_file), conf)
	

# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	if len(sys.argv) < 2:
		raise IndexError("Argument for config file is missing!")
	
	conf = sir.read_config(sys.argv[1])

	if len(sys.argv) < 3:
		raise IndexError("Argument for x is missing!")
	
	if conf['mode'] == "MC":
		if len(sys.argv) < 4:
			raise IndexError("Argument for savepath is missing!")
		savepath = sys.argv[3]
	else:
		if len(sys.argv) < 4:
			raise IndexError("Argument for y is missing!")
		if len(sys.argv) < 5:
			raise IndexError("Argument for savepath is missing!")
		
		savepath = sys.argv[4]
	
	if conf['mode'] == "MC":
		extract_profile_model_MC(conf, int(sys.argv[2]), savepath)
	elif conf['mode'] == "1C":
		extract_profile_model_1C(conf, int(sys.argv[2]), int(sys.argv[3]), savepath)
	elif conf['mode'] == "2C":
		extract_profile_model_2C(conf, int(sys.argv[2]), int(sys.argv[3]), savepath)



