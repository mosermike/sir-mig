# SIR MG

## Name
 SIR Multiple Guesses

## Version
main

## Description
The code computes physical parameters such as magnetic field, azimuth, inclination by using the SIR code. This code is used to perform inversions with creating random initial guesses. There are three different modes available:
 - MC : Monte Carlo Simulation
 - 1C : 1 component inversion
 - 2C : 2 components inversion

The mode is selected with the option "mode : XX" in the config file

Each python script has a help page which can be printed by using '-h' as an argument.

There is no introduction into SIR. For more information, look into the manual of the code.

**If you use this code, please cite these papers:**

## Installation
Several packages are necessary to use the sripts:
- numpy
- matplotlib
- scipy
- astropy

I used anaconda for executing the python scripts with the command python.

To install the SIR code, use

```
make fc=gfortran sir.x

```
in the src directory of the SIR code.

## Usage
1. To get started, copy an example directory to wherever you want. Adapt or create the config file with the existing python script create_config.py or by running
	```
		python [Path to here]/src/create_config.py
	```
   and the python script will ask you questions. One can also add the path, where the config file is saved, with executing the python script with '-path [Path]'. Do not use tabs in the config file, only spaces!

2. There is a python script 'main.py' which merges, normalises the data and corrects it for the spectral veil. Afterwards, the inversion is executed.

	```
		mpirun -np 64 python [Path to here]/src/main.py config.txt -dir [Path to the data directory of the fits files]
	```

   It executes for example the following commands for the modes 1C and 2C:

  1. Merging the data can be done as
	```
		python [path to here]/src/merge.py [config file] [directory of the fits, ending with /] [optional determination of the dataset number 001 (standard),002,003,004 for GRIS]
	```
	and in addition, some informations are saved and the wavelength numpy file is created
	

  2. Normalising the data with
	```
		python [path to here]/src/normalise.py [config file]
	```
  3. Correcting for the spectral veil (not for hinode):
	```
		python [path to here]/src/correction_spectral_veil.py [config file]
	```

  4. Performing the inversion:
	```
		mpirun -np [num. of processes] python [path to here]/src/inversion.py [config file]
	```
 3. To plot the results at specific wavelengths and log taus, use the script
	```
		python [path to here]/src/plot/result.py [config.py]
	```
	Use it with -h to see all the options.

 4. To plot the Stokes Profiles and the Model in one pixel, one can use the visualizer
	```
		python [path to here]/src/plot/visualizer.py [config file] [wavelength position]
	```
	use the flag '-h' for all the options.

## Example

In the directory `example`, several files are put. An examplary config file for gris and hinode for each mode. Then there are also other things needed to perform the inversion. All other files are created by using the scripts.

## Directory structure of src

```bash
.
├── create_config.py                  - Create config file as expected from the code
├── definitions.py                    - Global definitions
├── main.py                           - Main script to be executed
├── mml.mplstyle                      - Matplotlib stylesheet
├── model_1C.py                       - Class for the model of 1 Component
├── model_2C.py                       - Class for the model of 2 Components (includes the filling factor)
├── obs.py                            - Scripts related to observations
├── profile_stk.py                    - Class profile 
├── sir.py                            - Script containing import functions related to SIR and the code (reading config, writing control files, etc.)
├── _1C                               - Code for the 1 component inversion
│   ├── create_random_guess.py        - Creates random guesses for the inversion
│   └── inversion.py                  - Performs the inversion
├── _2C                               - Code for the 2 component inversion
│   ├── create_random_guess.py        - Creates random guesses for the inversion
│   └── inversion.py                  - Performs the inversion
├── _MC                               - Code for the Monte-Carlo-Simulation
│   ├── add_noise.py                  - Adding noise to the created, random profiles
│   ├── create_models.py              - Create models with 2 or 3 nodes and T based on HSRA.
│   ├── create_random_guess.py        - Creates random guess models
│   ├── inversion.py                  - Performs the MC inversion
│   └── synthesis.py                  - Performs the MC synthesis
├── plots                             - Directory with several plots for the different modes
│   ├── _1C
│   │   ├── inversion.py              - Plots one single inversion and saves it by using the config file
│   │   ├── inversion_single.py       - Plots one single inversion and saves it
│   │   └── result.py                 - Plots the results of the inversion
│   ├── _2C
│   │   ├── inversion.py              - Plots one single inversion and saves it by using the config file
│   │   ├── inversion_single.py       - Plots one single inversion and saves it
│   │   └── result.py                 - Plots the results of the inversion
│   ├── _MC
│   │   ├── analysis_compare_chi2.py  - Compares the chi^2 value for different MC simulations. Needs to be revised.
│   │   ├── analysis_multiple.py	  - Does the same as the script below but plots the results of multiple different MC simulations.
│   │   ├── analysis.py               - Analysis the MC simulation by plotting the standard deviations and printing out the values at some log tau values.
│   │   ├── chi2.py                   - Plots the chi2 values of multiple MC simulations 
│   │   ├── inversion.py              - Plots the Stokes profiles and the fit plus the physical parameters of one model.
│   │   ├── inversion_2.py            - Does the same as the script before but for two inversions.
│   │   ├── inversion_single.py       - Plots a single inversion result without any config file for testing purposes.
│   │   ├── synthesis_blend.py        - Blends two profiles with specified values alpha. It is also possible to add noise to the profiles.
│   │   └── synthesis.py              - Plots the results of the synthesis for multiple profiles and models.
|   ├── Ic_test.py                    - Plots the intensity of the obs. and the fit in a diagramm to test the inversion
│   ├── initial_vs_result.py          - Plots the initial starting point vs the resulting point for different inversions
|   ├── gris_sketch.py                - Sketches data of gris to get an overview
│   ├── one_model.py                  - Plots one model
│   ├── response.py                   - Plots the response function of a model created with SIR
│   └── visualizer.py                 - Visualises the Stokes Profiles and the fits as well as the model for selectable pixels.
├── preprocess                        - Scripts for data preprocessing
│   ├── correction_spectral_veil.py   - Correcting the spectral veil for ground based telescopes by using FTS data
│   ├── merge.py                      - Merging hinode and gris data to a data cube
│   └── normalise.py                  - Normalises Hinode and GRIS data based on the defined wavelength in definitions.py
└── tools                             - Useful scripts
    ├── change_config_path.py         - Changes the path of a config file by providing the new path in the command line
    ├── change_model.py               - Changes a parameter of a model. Verify that it worked as wanted.
    ├── clean.py                      - Cleans the inversion files
    ├── convert_MC.py                 - A script to convert old MC files to the new version (separated profiles saved for different lines)
    ├── distclean.py                  - Removes all created files
    └── extract_profile_model.py      - Extract the data from the npy files and stores it in a directory to run a single inversion again


```

## Note
- Example models can be found in the directory 'Models'
- An example structure which is expected can be found in Example
- For the input data, the code assumes the unnormalised data. Normalisation and spectral veil correction can then be done and a prefix is added automatically (defined in definitions.py) and the corrected data is used for the inversion! 
- If you do not want to correct for the spectral veil but you want to use a PSF and let it be created by the code, do the following:
  1. Create a numpy array with one entry. This entry will be taken as the width of a Gaussian in mA.
  2. Save this numpy array as a npy file with the same name as the parameter 'veil_parameters' in the definitions python script in 'src' (if it is not changed, it should be 'veil.npy').
  3. make sure that there is no file in the directory with the same name as provided in the config file

## Support
For any questions, write a mail to [mosermike@protonmail.ch](mailto:mosermike@protonmail.ch)


## Authors and acknowledgment
I thank Dr. Juan Manuel Borrero and Dr. Ivan Milic for helpful insights into the topic. I also thank Dr. Juan Manuel Borrero for giving me the inversion code.


## Project status
Ongoing


## License
This software is distributed under GPL 2.0 license. You can use/modifiy/distribute as you wish as long
as it is under the GPL 2.0 license.
