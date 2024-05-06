"""
Here are some global definitions stored. This are used during different processes. They are needed to identify specific outputs.
You can change them as you want. They should not have any real impact but different naming. Only the parameter ll_lit has an impact.
Ideally, they are not changed. Keep in mind that the part in data stuff is used for plotting results.

I DO NOT RECOMMEND FOR CONSISTENCY TO CHANGE THOSE PARAMETERS WITH THE COMMENT 'KEEP' WHEN SIMULATIONS WERE ALREADY PERFORMED.
THIS COULD RESULT INTO ERRORS IF FILES ARE READ AT A LATER POINT.

PARAMETERS WITH THE COMMENT 'CLEAN' WILL NOT BE DELETED IF CHANGED AFTER A SIMULATION IS RUN AND 'distclean.py' OR 'clean.py' IS
EXECUTED

"""

##############
# Data stuff #
##############
end_stokes = '_stokes.bin' # Ending inversion result stokes  KEEP
end_models = '_models.bin' # Ending inversion result models  KEEP
end_models1 = '_models_1.bin'  # Ending inversion result models 1 for 2C KEEP
end_models2 = '_models_2.bin'  # Ending inversion result models 2 for 2C KEEP 
end_norm = '_norm.bin'  # Normalised ending in file KEEP
end_errors = '_errors.bin'  # Ending inversion result error KEEP
end_errors1 = '_errors_1.bin'  # Ending inversion result error 1 for 2C KEEP
end_errors2 = '_errors_2.bin'  # Ending inversion result error 2 for 2C KEEP
veil_parameters = 'veil.npy'  # numpy file with sigma and nu from the spectral veil correction KEEP
header_infos = 'infos.txt'  # Information about the observations from the Header # KEEP
filling_factor = 'fill'		# Generic start of the npy files with the filling factors (end_modelsX is added)

#####################
# Synthesis Stuff   #
#####################
model_syn = 'syn.mod' # Name of the synthesis model
model_inv = 'guess.mod' # Name of the inversion guess model
syn_trol_file ='syn.trol' # inversion control file CLEAN

#####################
# Inversion stuff	#
#####################
Grid = 'Grid.grid'  # Grid file KEEP
profile_obs = 'profile.per'  # Name of the observed profile in the inversion
inv_trol_file = 'inv.trol'  # inversion control file CLEAN
guess = "guess.mod"  # Guess model as in control file for 1C and MC
guess1 = 'guess1.mod'  # Guess Model 1 for control file and inversion for 2C
guess2 = 'guess2.mod'  # Guess Model 2 for control file and inversion for 2C
best_guess = 'best_guess.mod'  # Name of the best guess model (best model is saved under this name) KEEP
best_guess1 = 'best_guess1.mod'  # Name of the best guess model (best model is saved under this name) for 2C KEEP
best_guess2 = 'best_guess2.mod'  # Name of the best guess model (best model is saved under this name) for 2C KEEP
task_start = ".task_"	# Start of the task folders
profile = 'profile.per' # Name of the observed profile in the synthesis/inversion
model = 'model_'  # Generic name for creating guesses for Model 1 for 1C
model1 = 'model1_'  # Generic name for creating guesses for Model 1 for 2C
model2 = 'model2_'  # Generic name for creating guesses for Model 2 for 2C
profile_obs = profile  # Name of the observed profile in the inversion
header_2C = "   0.10000000      0.50000	   0.00000" # Header used for writing a model 2C
fill_2C = 0.5 # Filling factor for the first model when mode 2C is used

#############################
# Spectral Veil Correction	#
#############################
# Literature wavelength in air
# Add here for another instrument
ll_lit = {
			"GRIS": 15648.514  # From A new Multiplet table for Fe I
		}


#####################
# Normalisation		#
#####################
# Dict. describing the range where the continuum is expected for normalisation.
# If you want another instrument to be implemented, add it here.
# The script merge might need to be adapted manually.
ll_lit_norm = {
			"GRIS": [15640.5, 15643],
			"Hinode": [6301, 6301.072]
			}

################
# For plotting #
################
origin = "lower"
# Add here another instrument if wanted for plotting in a different wavelength value (in rel. values)
# This value is subtracted from the wavelength
ll_relative = {
				"GRIS": 15600,
				"Hinode": 6300
}

#################
# Miscellaneous	#
#################
# All defined instruments, add your instrument here if you want to add another one
# Merge needs probably adaption. Especially where I noted ADAPT
instruments = ["GRIS", "Hinode"]
plt_lib = "" # Name for the matplotlib style sheet to plot if not the mml.mplstyle should be used, use the absolute path


