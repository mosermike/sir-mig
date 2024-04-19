# SIR-MIG-MC

## Name
 SIR Monte Carlo Simulation with Multiple Initial Guesses

## Version
main

## Description
The code computes physical parameters such as magnetic field, azimuth, inclination by using the SIR code. This code is intended to first create random models, then perform a synthesis. Afterwards, noise is added to all profiles. In the end, the noisy profiles are inverted either by a specific guess model or by random guess models. The intentip=on is to perform a Monte-Carlo simulation.

Each python script has a help page which can be printed by using '-h' as an argument.

There is no introduction into SIR. For more information, look into the manual of the code.

**If you use this code, please cite these papers:**

## Installation
Several packages are necessary to use the sripts:
- numpy
- matplotlib
- scipy

I used anaconda for executing the python scripts with the command python.

To install the SIR code, use

```
make fc=gfortran sir.x

```
in the src directory of the SIR code.

## Usage
1. To get started, copy an example directory to wherever you want. Adapt or create the config file with the existing python script create_config.py by running
	```
		python [path to sir-mc-model]/src/tools/create_config.py
	```
and the python script will ask you questions. One can also add the path, where the config file is saved, with executing the python script with '-path [Path]'. Do not use tabs in the config file, only spaces!

2. There is a main script which performs all 4 parts of the Monte Carlo simulation with SIR and different synthesised models.
   The four parts are the following:
	- Create Models
	- Synthesise the models
	- Add Noise to the Profiles
	- Perform the inversion

To perform the MC from scratch, perform:
	```
		mpirun -np [Number of parallel processes] python src/main.py [config file]
	```

It executes following commands:

 1. Creating Models:
	```
		python src/create_models.py [config file]
	```
	with 2 or 3 nodes. This is determined by the config file.

 2. Synthesising the Models:
	```
		mpirun -np [# processes] python src/synthesis.py [config file]
	```
 3. Adding noise to the profiles:
	```
		python src/add_noise.py [config file]
	```

 4. Performing the inversion:
	```
		mpirun -np [parallel jobs] python src/inversion.py [config file] 
	```
	If this process takes forever, this is probably caused by the i/o ability of the hard drive.


To analyse the MC simulation (computing the deviation from the inversion to the simulated model), one can use
```
		python src/plot/analysis.py [config file]
```
Use the command with -h to see a help page for the script.

## Example

In the directory `example`, several files are put. A model which the synthesis is based on, a grid and a line file for the SIR inversion, and two control files one for the synthesis and one for the inversion for data as given from GRIS and for data as given from Hinode. For example, a fast MC simulation based on the example can be done by doing the following:
1. Execute
```
	python [path to  sir-mc-model]/src/change_config_path.py config.txt $PWD
```
in the directory of the config.txt file. 
2. Start the MC simulation with

	```
		mpirun -np [# processes] python [path to sir-mc-model]/src/main.py [config file]
	```
	and wait for the result.

## Directory structure
A typical directory with a finished simulation looks like that:
```bash
.
├── model.mod - Basis model
├── config.txt - Config file
├── Grid.grid - Grid determining which lines 
├── inversion*.npy - Inversion results (Models, Errors, Stokes)
├── chi2.npy - Chi2 map
├── inv.trol - Control file for the inversion
├── Lines - Line files
├── noise_profiles_\*.npy - Noise profiles as npy file (\* stands for the line number)
├── sir.x - Executable sir file
├── syn_profiles_\*.npy - Synthesis profiles (\* stands for the line number)
├── syn.trol - Control file for the synthesis
└── THEVENIN - Abundance file from the SIR code

```

## Note
- Example models can be found in the directory 'Models'
- An example structure which is expected can be found in Example
- The path parameter in the config file determines where the directories are expected like syn, noise and inversion as well as the SIR files!
- The first two entries in the 'create_limX' settings in the config file are used for the node deeper in the atmosphere (e.g. log tau 1) and the second one higher up (e.g. log tau -3).

## Python Scripts
Here, the functions of the scripts in the folder src are explained:
```bash
.
├── add_noise.py					- Adding noise to the created, random profiles
├── create_models.py				- Create models with 2 or 3 nodes and T based on HSRA.
├── create_random_guess.py			- Creates random guess models
├── inversion.py					- Performs the MC inversion with mpirun
├── lib							- Scripts used in (nearly) all scripts
│   ├── definitions.py				- Predefined names and values
│   ├── model.py					- File for the class Model and Error
│   └── sir.py						- Library for different functions related to SIR
├── main.py						- Script which performs every step of the MC simulation
├── plot							- Scripts for plotting stuff
│   ├── analysis_compare_chi2.py		- Compares the chi^2 value for different MC simulations. Needs to be revised.
│   ├── analysis_multiple.py			- Does the same as the script below but plots the results of multiple different MC simulations.
│   ├── analysis.py					- Analysis the MC simulation by plotting the standard deviations and printing out the values at some log tau values.
│   ├── chi2.py					- Plots the chi2 values of multiple MC simulations
│   ├── initial_vs_result.py			- Plots the initial start position and the result of one inversion at a given log tau.
│   ├── inversion.py				- Plots the Stokes profiles and the fit plus the physical parameters of one model.
│   ├── inversion_2.py				- Does the same as the script before but for two inversions.
│   ├── inversion_single.py			- Plots a single inversion result without any config file for testing purposes.
│   ├── mml.mplstyle				- Library for plotting stuff
│   ├── one_model.py				- Plots the physical parameters of one model.
│   ├── response.py					- Plots the response functions created with SIR
│   ├── synthesis_blend.py			- Blends two profiles with specified values alpha. It is also possible to add noise to the profiles.
│   └── synthesis.py				- Plots the results of the synthesis for multiple profiles and models.
├── synthesis.py					- Performs the MC synthesis with mpirun
└── tools							- Scripts for different things to help
    ├── change_config_path.py			- Changes the config path in a config file
    ├── change_model.py				- Changes a model to defined values by multiplying or adding a factor.
    ├── clean.py					- Cleans up the inversion.
    ├── create_config.py				- Creating the config file by asking questions.
    ├── distclean.py				- Cleans up all created files but pictures and config file.
    └── extract_profile_model.py		- Extracts the profiles and models from the npy files for a specific model
```

## Support
For any questions, write a mail to [mikemoser1997@gmail.com](mailto:mikemoser1997@gmail.com)


## Authors and acknowledgment
I thank Dr. Juan Manuel Borrero and Dr. Ivan Milic for helpful insights into the topic. I also thank Dr. Juan Manuel Borrero for giving me the inversion code.


## Project status
Ongoing


## License
This software is distributed under GPL 2.0 license. You can use/modifiy/distribute as you wish as long
as it is under the GPL 2.0 license.
