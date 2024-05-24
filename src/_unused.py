######################################################################

def write_profile(filename, profiles, pos):
	"""
	Write a profile for a specific model number to a file

	Parameters
	----------
	filename : string
		Name of the saved file
	profiles : list
		List containing all the profiles
	atoms : list
		List containing the number of the line from the Line file
	pos : int
		Position which model is saved

	Returns
	-------
	None
	"""

	f = open(filename, 'w')
	#for a in range(len(atoms)):
	for m in range(profiles.shape[2]):
			f.write(f" {int(profiles[pos,0,m]):>2} {profiles[pos,1,m]:>10.4f} {profiles[pos,2,m]:>14.7E} {profiles[pos,3,m]:>14.7E} {profiles[pos,4,m]:>14.7E} {profiles[pos,5,m]:>14.7E}\n")
	f.close()
