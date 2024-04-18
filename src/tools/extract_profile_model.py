"""
Extracts the profile file and the mdoel file from a data cube
"""
import numpy as np
import sys, os, shutil
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import obs, sir, definitions as d
from os.path import exists
import profile_stk as p
import model_1C as m
import model_2C as m2

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
	obs1 = p.Profile(os.path.join(path, conf['cube_inv']),os.path.join(path, conf['waves']))
	obs1.cut_to_wave(sir.angstrom_to_pixel(obs1.wave, conf["range_wave"])) # Cut the values to data used in the inversion

	inv = p.Profile(os.path.join(path, conf['inv_out']) + d.end_stokes,os.path.join(path, conf['inv_out'])+d.end_wave)
	if(inv.nw == obs1.nw):
		inv.cut_to_wave(sir.angstrom_to_pixel(inv.wave, conf["range_wave"])) # Cut the values to data used in the inversion
	mod = m.Model(filename=os.path.join(path,conf['inv_out'] + d.end_models))
	guess = m.Model(filename=os.path.join(path,d.best_guess.replace(".mod",".npy")))
	err = m.Error(filename=os.path.join(path,conf['inv_out'] + d.end_errors))

	# Change to the reduced Map
	x = x - Map[0]
	y = y - Map[2]
	
	# Save data in formats for SIR
	obs1.write(savepath + "profile" + add + ".per", x + Map[0], y + Map[2], d.Grid)
	inv.write(savepath + "profile_result" + add + ".per", x, y, d.Grid)
	mod.write(savepath + "model_result" + add + ".mod", d.header, x, y)
	guess.write(savepath + conf["model"].replace(".mod","") + add + ".mod", d.header, x, y)
	err.write(savepath + "model_result" + add + ".err", d.header, x, y)
	
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
	
	Parameter
	---------
	config : dict
		Dictionary containing all the information from the config
	num : int
		Model number

	Return
	------
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
	# Write the data from the cube into a profile file for SIR
	atoms = [i.split(',') for i in conf['atoms']]
	noise_profiles = []
	syn_profiles = []
	inv_profiles = []
	noise = p.Profile_MC(os.path.join(path,conf['noise_out']))
	syn = p.Profile_MC(os.path.join(path,conf['syn_out']))
	inv = p.Profile_MC(os.path.join(path,f"{conf['inv_out']}{d.end_stokes}"))

	noise.write(os.path.join(savepath,'noise.per' + add), num-1)
	syn.write(os.path.join(savepath,'syn.per' + add), num-1)
	inv.write(os.path.join(savepath,'inv.per' + add), num-1)
	

	# Now with the class model/error
	mod = m.Model(os.path.join(path,conf['inv_out'] + d.inv_models))
	gue = m.Model(os.path.join(path,d.best_guess.replace('.mod','.npy')))
	syn = m.Model(os.path.join(path,conf['model_out']))
	err = m.Error(os.path.join(path, os.path.join(path,conf['inv_out'] + d.inv_errors)))
	
	syn.write(os.path.join(savepath,'syn.mod' + add), d.header,num-1)
	mod.write(os.path.join(savepath,'res.mod' + add), d.header,num-1)
	gue.write(os.path.join(savepath,'guess.mod' + add), d.header,num-1)
	err.write(os.path.join(savepath,'res.err' + add), d.header,num-1)

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

	# Load wavelenghts and map of (x,y) used in the inversion
	Map = conf['map']
	waves = np.load(os.path.join(path, conf['waves']))
	
	# Load data
	obs1 = obs.load_data(conf, filename=conf['cube_inv'])	
	inv = np.load(os.path.join(path,conf['inv_out'] + d.end_stokes))
	guess1 = np.load(os.path.join(path,d.best_guess1.replace(".mod",".npy")))
	guess2 = np.load(os.path.join(path,d.best_guess2.replace(".mod",".npy")))
	mod1 = np.load(os.path.join(path,conf['inv_out'] + d.end_models1))
	mod2 = np.load(os.path.join(path,conf['inv_out'] + d.end_models2))
	err1 = np.load(os.path.join(path,conf['inv_out'] + d.end_errors1))
	err2 = np.load(os.path.join(path,conf['inv_out'] + d.end_errors2))
	

	# Cut observation
	obs1 = obs1[Map[0]:Map[1]+1,Map[2]:Map[3]+1, :, :]

	# Change to the reduced Map
	x = x - Map[0]
	y = y - Map[2]
	
	# Save data in formats for SIR
	obs.write_profile(savepath + "profile" + add + ".per", obs1, conf, x, y)
	obs.write_profile(savepath + "profile_result" + add + ".per", inv, conf, x, y)
	sir.write_model(savepath + "guess1" + add + ".mod", d.header,
						 guess1[x,y,0],guess1[x,y,1],guess1[x,y,2],guess1[x,y,3],
						 guess1[x,y,4],guess1[x,y,5],guess1[x,y,6],guess1[x,y,7],
						 guess1[x,y,8],guess1[x,y,9],guess1[x,y,10]
						)
	sir.write_model(savepath + "guess2" + add + ".mod", d.header,
						 guess2[x,y,0],guess2[x,y,1],guess2[x,y,2],guess2[x,y,3],
						 guess2[x,y,4],guess2[x,y,5],guess2[x,y,6],guess2[x,y,7],
						 guess2[x,y,8],guess2[x,y,9],guess2[x,y,10]
						)
	sir.write_model(savepath + "model_result1" + add + ".mod", d.header,
						 mod1[x,y,0],mod1[x,y,1],mod1[x,y,2],mod1[x,y,3],
						 mod1[x,y,4],mod1[x,y,5],mod1[x,y,6],mod1[x,y,7],
						 mod1[x,y,8],mod1[x,y,9],mod1[x,y,10]
						)
	sir.write_model(savepath + "model_result2" + add + ".mod", d.header,
						 mod2[x,y,0],mod2[x,y,1],mod2[x,y,2],mod2[x,y,3],
						 mod2[x,y,4],mod2[x,y,5],mod2[x,y,6],mod2[x,y,7],
						 mod2[x,y,8],mod2[x,y,9],mod2[x,y,10]
						)
	sir.write_model(savepath + "model_result1" + add + ".err", d.header,
						 err1[x,y,0],err1[x,y,1],err1[x,y,2],err1[x,y,3],
						 err1[x,y,4],err1[x,y,5],err1[x,y,6],err1[x,y,7],
						 err1[x,y,8],err1[x,y,9],err1[x,y,10]
						)

	sir.write_model(savepath + "model_result2" + add + ".err", d.header,
						 err2[x,y,0],err2[x,y,1],err2[x,y,2],err2[x,y,3],
						 err2[x,y,4],err2[x,y,5],err2[x,y,6],err2[x,y,7],
						 err2[x,y,8],err2[x,y,9],err2[x,y,10]
						)

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



