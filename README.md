# SIR MG
	
## Name
 SIR Multiple Guesses

## Version
1.0

## Description
The code computes physical parameters such as magnetic field, azimuth, inclination by using the SIR code. This code is used to perform inversions with creating random initial guesses. There are three different modes available:
 - MC : Monte Carlo Simulation
 - SY : Synthesis
 - 1C : 1 component inversion
 - 2C : 2 components inversion

The mode is selected with the option "mode : XX" in the config file

Each python script has a help page which can be printed by using '-h' as an argument.

There is no introduction into SIR. For more information, look into the manual of the code.

**There is a wiki available at [sir-mig.moser.mywire.org](https://sir-mig.moser.mywire.org) **


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
├── chi2_stk.py                       - Class for computing chi2
├── create_config.py                  - Create config file as expected from the code
├── create_random_guess.py            - Creates random guesses for the inversion
├── definitions.py                    - Global definitions
├── inversion.py                      - Performs the inversion
├── main.py                           - Main script to be executed
├── mml.mplstyle                      - Matplotlib Style Sheet
├── model_atm.py                      - Class for the model
├── preprocess.py                     - Scripts for data preprocessing (merge, normalise, spectral veil correction)
├── profile_stk.py                    - Class profile 
├── simulation.py                     - Performs the simulation which creates models and stokes profiles.
├──sir.py                             - Script containing import functions related to SIR and the code (reading config, writing control files, etc.)
└──synthesis.py                       - Perform synthesis for mode "SY"

```

## Note
- Example models can be found in the directory 'Models'
- An example structure which is expected can be found in Example
- For the input data, the code assumes the unnormalised data. Normalisation and spectral veil correction can then be done and a prefix is added automatically (defined in definitions.py) and the corrected data is used for the inversion! 
- If you do not want to correct for the spectral veil but you want to use a PSF and let it be created by the code, do the following:
  1. Create a numpy array with one entry. This entry will be taken as the width of a Gaussian in mA.
  2. Save this numpy array as a npy file with the same name as the parameter 'veil_parameters' in the definitions python script in 'src' (if it is not changed, it should be 'veil.npy').
  3. make sure that there is no file in the directory with the same name as provided in the config file


## FTS data
	FTS data for the spectral veil correction can be found [here](https://aasfiles.blob.core.windows.net/files/archive/cdrom/volume7/doc/files7.htm#1996ApJS..106..165W)


## Support
For any questions, write a mail to [mosermike@protonmail.ch](mailto:mosermike@protonmail.ch)


## Authors and acknowledgment
I thank Dr. Juan Manuel Borrero and Dr. Ivan Milic for helpful insights into the topic. I also thank Dr. Juan Manuel Borrero for giving me the inversion code.


## Project status
Finished, only updated are made.


## License
This software is distributed under GPL 2.0 license. You can use/modifiy/distribute as you wish as long
as it is under the GPL 2.0 license.
