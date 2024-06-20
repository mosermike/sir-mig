"""
Here are some global definitions stored. This are used during different processes. They are needed to identify specific outputs.
You can change them as you want. They should not have any real impact but different naming. Only the parameter ll_lit has an impact.
Ideally, they are not changed. Keep in mind that the part in data stuff is used for plotting results.

I DO NOT RECOMMEND FOR CONSISTENCY TO CHANGE THOSE PARAMETERS WITH THE COMMENT 'KEEP' WHEN SIMULATIONS WERE ALREADY PERFORMED.
THIS COULD RESULT INTO ERRORS IF FILES ARE READ AT A LATER POINT.

PARAMETERS WITH THE COMMENT 'CLEAN' WILL NOT BE DELETED IF CHANGED AFTER A SIMULATION IS RUN AND 'distclean.py' OR 'clean.py' IS
EXECUTED

"""
import numpy as np

##############
# Data stuff #
##############
end_stokes = '_stokes.bin' # Ending inversion result stokes  KEEP
end_models = '_models.bin' # Ending inversion result models  KEEP
end_models1 = '_models1.bin'  # Ending inversion result models 1 for 2C KEEP
end_models2 = '_models2.bin'  # Ending inversion result models 2 for 2C KEEP 
end_errors = '_errors.bin'  # Ending inversion result error KEEP
end_errors1 = '_errors1.bin'  # Ending inversion result error 1 for 2C KEEP
end_errors2 = '_errors2.bin'  # Ending inversion result error 2 for 2C KEEP
veil_parameters = 'veil.npy'  # numpy file with sigma and nu from the spectral veil correction KEEP
header_infos = 'infos.txt'  # Information about the observations from the Header # KEEP
best_guess_file = '_best_guess.bin'  # Name of the best guess model (best model is saved under this name) KEEP
best_guess1_file = '_best_guess1.bin'  # Name of the best guess model (best model is saved under this name) for 2C KEEP
best_guess2_file = '_best_guess2.bin'  # Name of the best guess model (best model is saved under this name) for 2C KEEP
cube = 'data.bin' # Data cube from the fits files
cube_norm = 'data_norm.bin' # Normalised data cube

#####################
# Synthesis Stuff   #
#####################
model_syn = 'syn.mod' # Name of the synthesis model
model_inv = 'guess.mod' # Name of the inversion guess model
syn_trol_file ='syn.trol' # inversion control file CLEAN

#####################
# Inversion stuff	#
#####################
Grid = 'Grid.grid'  # Grid file
profile_obs = 'profile.per'  # Name of the observed profile in the inversion
inv_trol_file = 'inv.trol'  # inversion control file CLEAN
guess = "guess.mod"  # Guess model as in control file for 1C and MC
guess1 = 'guess1.mod'  # Guess Model 1 for control file and inversion for 2C
guess2 = 'guess2.mod'  # Guess Model 2 for control file and inversion for 2C
best_guess = 'best_guess.mod'  # Name of the best guess model  (do not choose best.mod)
best_guess1 = 'best_guess1.mod'  # Name of the best guess model for 2C (do not choose best1.mod)
best_guess2 = 'best_guess2.mod'  # Name of the best guess model for 2C  (do not choose best2.mod)
task_start = ".task_"	# Start of the task folders
profile = 'profile.per' # Name of the observed profile in the synthesis/inversion
model = 'model_'  # Generic name for creating guesses for Model 1 for 1C
model1 = 'model1_'  # Generic name for creating guesses for Model 1 for 2C
model2 = 'model2_'  # Generic name for creating guesses for Model 2 for 2C
profile_obs = profile  # Name of the observed profile in the inversion
header_2C = "   0.10000000      0.50000	   0.00000" # Header used for writing a model 2C
fill_2C = 0.5 # Filling factor for the first model when mode 2C is used
psf = ".psf.dat" # PSF file for the inversion

#############################
# Spectral Veil Correction	#
#############################
# Literature wavelength in air
# Add here for another instrument
ll_lit = {
			"GRIS": 15648.514  # From A new Multiplet table for Fe I
		}

sigma_range = (10,150) # in mA
nu_range    = (0,31) # in %


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
# All defined instruments, add your instrument here if you want to add another one. Add it at the end and not between the strings or at the start
# Merge needs probably adaption. Especially where I noted ADAPT
instruments = ["GRIS", "Hinode"]
plt_lib = "" # Name for the matplotlib style sheet to plot if not the mml.mplstyle should be used, use the absolute path


###############################
# Create Temperature Profiles #
###############################
###########################################################################################
# The two models hsra and cool11 (Collados M., Martínez Pillet V., Ruiz Cobo B., 		  #
# Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 (Umbral model for a big spot) are  #
# considered. The minimum is the cool11 model and then I add a factor of the HSRA model.  #
# The structure of this factor is the following:									      #
# The factors are chosen because of the following thoughts/points:					      #
# - The temperature starts at the model cool11.									          #
# - I create a random factor between 0 and 1.									          #
# - 0 = cool Model																		  #
# - 1 = HSRA Model														  	              #
# - I implemented some dependency on the magnetic field: If the magnetic field is strong, #
#   the range for the factor is smaller											          #
###########################################################################################
log_taus = np.array([ 1.4,  1.3,  1.2,  1.1,  1. ,  0.9,  0.8,  0.7,  0.6,  0.5,  0.4,
						   0.3,  0.2,  0.1,  0. , -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7,
						  -0.8, -0.9, -1. , -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8,
						  -1.9, -2. , -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9,
						  -3. , -3.1, -3.2, -3.3, -3.4, -3.5, -3.6, -3.7, -3.8, -3.9, -4. ])

# HSRA Model for the upper boundary
upper_T = np.array([9560. , 9390. , 9220. , 9050. , 8880. , 8710. , 8520. , 8290. ,
					  8030. , 7750. , 7440. , 7140.9, 6860. , 6610. , 6390. , 6200. ,
					  6035. , 5890. , 5765. , 5650. , 5540. , 5430. , 5330. , 5240. ,
					  5160. , 5080. , 5010. , 4950. , 4895. , 4840. , 4790. , 4750. ,
					  4720. , 4690. , 4660. , 4630. , 4600. , 4575. , 4550. , 4525. ,
					  4490. , 4460. , 4430. , 4405. , 4380. , 4355. , 4330. , 4305. ,
				       4280. , 4250. , 4225. , 4205. , 4190. , 4175. , 4170. ])

# cool penumbra model for lower boundary
# (Collados M., Martínez Pillet V., Ruiz Cobo B.,Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 (Umbral model for a big spot)
lower_T = np.array([6780.3, 6536. , 6291.9, 6048.5, 5806.5, 5569.5, 5340.7, 5117.3,
					  4902.9, 4700.4, 4513.9, 4342.3, 4188.8, 4053.4, 3940.5, 3854. ,
					  3785. , 3726.8, 3676.7, 3633.6, 3597.9, 3564.7, 3534.9, 3511.6,
					  3498. , 3489.4, 3482.2, 3475.6, 3468.9, 3461.6, 3453.6, 3445.2,
					  3436.4, 3427. , 3417.1, 3406.5, 3395.3, 3383.4, 3370.8, 3357.9,
					  3345.1, 3332.4, 3319.2, 3305.5, 3291.1, 3276. , 3260.1, 3243.5,
					  3225.9, 3207.5, 3188.5, 3170.5, 3155.7, 3142.8, 3129.7]
					)

# Multiplicative random factor (the cool model is multiplied with a number between 1-f and 1+f)
multiplicative_T = 0.05

