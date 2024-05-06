"""
Misc. Functions
"""

def initial(mode):
	"""
	Initial print outs and preparation

	Parameters
	----------
	mode : string
		Mode which is used
	
	Returns
	-------
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

def x_y_add_zeros(x, y):
	"""
	Adds zeros so that the returning strings have 4 letters

	Parameters
	----------
	x : float
		x position
	y : float
		y position

	Returns
	-------
	out : str
		x as a string of 4 letters
	out : str
		y as a string of 4 letters

	"""
	if x < 10:
		x_str = "000" + str(x)
	elif x < 100:
		x_str = "00" + str(x)
	elif x < 1000:
		x_str = "0" + str(x)
	else:
		x_str = str(x)
	if y < 10:
		y_str = "000" + str(y)
	elif y < 100:
		y_str = "00" + str(y)
	elif y < 1000:
		y_str = "0" + str(y)
	else:
		y_str = str(y)

	return x_str, y_str

def create_task_folder_list(arg):
	"""
	Creates a list which folders should be created and executed. This is done so
	that the inversion itself can be executed linearly to make use of all cores.

	Parameters
	----------
	arg : numpy array or int
		1x4 array containing the limits in x and y of the data or number of 1D models

	Returns
	-------
	out : dict
		Dictionary with all the names of the task folders, x and y position
	"""
	import numpy as np
	import definitions as d
	# Create arrays
	tasks = []
	xs = []
	ys = []
	
	if isinstance(arg,int):
		y = 0
		# Determine task folder names
		for x in range(arg):
			x_str, y_str = x_y_add_zeros(x,0)

			tasks.append(d.task_start + x_str + "_" + y_str)
			xs.append(x)
			ys.append(y)
		Dict = {
				'folders' : tasks,
				'x' : np.array(xs),
				'y' : np.array(ys),
		
		}
	else:
		Map = arg
		# Determine task folder names
		for x in range(Map[0], Map[1]+1):
			for y in range(Map[2], Map[3]+1):
				x_str, y_str = x_y_add_zeros(x, y)

				tasks.append(d.task_start + x_str + "_" + y_str)
				xs.append(x)
				ys.append(y)
		Dict = {
				'folders': tasks,
				'x': np.array(xs),
				'y': np.array(ys),
		
		}

	return Dict