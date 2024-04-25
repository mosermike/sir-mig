"""
Changes the path in config
"""

import sys
sys.path.append(sys.path[0] + "/..")
import sir


def help():
	"""
	Prints information about the script
	"""
	print("change_config_path - Changes the path in a config file while keeping the rest")
	print("Usage: python change_config_path.py [OPTION]")
	print()
	sir.option("[1. Pos.]","Input config file")
	sir.option("[2. Pos.]","New path")
	sys.exit()


def change_config_path(conf, path):
	"""
	Changes the path in the config file

	Parameter
	---------
	inp : string
		In- and output of the config file
	path : string
		New path in the config file
	"""

	inp = conf['filename']

	# Read config, change path and write it

	conf["path"] = path

	sir.write_config(inp, conf)


# Used if executed directly
if __name__ == "__main__":
	if "-h" in sys.argv:
		help()
	conf = sir.read_config(sys.argv[1], check=False, change_config=True)
	change_config_path(conf, sys.argv[2])
