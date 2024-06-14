=============
Configuration
=============

.. _config:

-------
General
-------

The first part of the configuration is very general and determines where and what is done.

+-----------+-------------------------------------------------------------------------------+
| Option    | Description                                                                   |
+===========+===============================================================================+
| mode      | Determine which part of the code is executed. Options are '1C', '2C' and 'MC' |
+-----------+-------------------------------------------------------------------------------+
| path      | Path location where everything is stored and will be saved                    |
+-----------+-------------------------------------------------------------------------------+

The next parts describe the options for the modes '1C' and '2C'. The reason for that is the big differences between these two modes and the mode 'MC'

-----------------------
Mode `1C` and Mode `2C`
-----------------------

Stuff from the data
===================

+------------+-------------------------------------------------------------------------------+
| Option     | Description                                                                   |
+============+===============================================================================+
| cube       | The binary file used for the inversion                                        |
+------------+-------------------------------------------------------------------------------+
| map        | Rectangle are to be used for the inversion (given as a list xmin,xmax,...)    |
+------------+-------------------------------------------------------------------------------+



Data Preprocessing
==================

+------------+-------------------------------------------------------------------------------+
| Option     | Description                                                                   |
+============+===============================================================================+
| preprocess | Determines whether data is preprocessed (1) or not (0)                        |
+------------+-------------------------------------------------------------------------------+
| instrument | Name of the instrument used (GRIS,Hinode or empty) to use predefined range    |
+------------+-------------------------------------------------------------------------------+
| quiet_sun  | Rectangle area used for normalisation and spectral veil correction            |
+------------+-------------------------------------------------------------------------------+
| fts_file   | File with FTS data to correct for spectral veil in the defined wavelengths    |
+------------+-------------------------------------------------------------------------------+
| shift_wave | Shift the wavelength array by this value to correct for shifts                |
+------------+-------------------------------------------------------------------------------+
| save_cube  | Save the preprocessed data (merged and normalised data cube) (1 or 0)         |
+------------+-------------------------------------------------------------------------------+

Notes:
- When you use the instrument `GRIS` the spectral veil and the normalisation is done in the range [15640.5, 15643] which can be changed in the script `definitions.py`. Consider adding another instrument called e.g. 'GRIS1' and add the wavelength range you want in case you do not consider data in the range 15600.
- When you use the instrument `Hinode`, the normalisation is done in the range [6301, 6301.072]. This can be changed in the script 'definitions.py
- When `fts_file` does not have an entry, no spectral veil correction is used.
- To select the fits files, use the flag '-dir [path to fits files]'  when you execute the main script.


Inversion configuration
=======================

+------------+----------------------------------------------------------------------------------+
| Option     | Description                                                                      |
+============+==================================================================================+
| model      | 1st Model Name used as base model for guesses (1C)                               |
+------------+----------------------------------------------------------------------------------+
| model1     | 1st Model Name used as base model 1 for guesses (2C)                             |
+------------+----------------------------------------------------------------------------------+
| model2     | 2nd Model Name used as base model 2 for guesses (2C)                             |
+------------+----------------------------------------------------------------------------------+
| range_wave | Range in the Grid (start in A, step in mA, number of wavelengths), `;` = newline |
+------------+----------------------------------------------------------------------------------+
| fill       | Filling factors for both models as `fill1,fill2` (if random_guess > 0) (2C)      |
+------------+----------------------------------------------------------------------------------+
| inv_out    | Prefix to written results from the inversion                                     |
+------------+----------------------------------------------------------------------------------+
| chi2       | Output of the chi2 values (npy)                                                  |
+------------+----------------------------------------------------------------------------------+
| line       | Lines File                                                                       |
+------------+----------------------------------------------------------------------------------+
| atoms      | Atoms used as a list where 1;1 = newline for the Grid File                       |
+------------+----------------------------------------------------------------------------------+
| guess      | Binary file used as initial guess, blank = base model                            |
+------------+----------------------------------------------------------------------------------+
| guess1     | Binary file used as 1st initial guess, blank = base model                        |
+------------+----------------------------------------------------------------------------------+
| guess2     | Binary file used as 2nd initial guess, blank = base model                        |
+------------+----------------------------------------------------------------------------------+
| psf        | Spectral psf file, 'gauss 1.0' = create Gaussian PSF with $\sigma = 1.0$         |
+------------+----------------------------------------------------------------------------------+


Notes:
- `Base Model` refers to the model where the data is taken which are not randomised
- For `guess(1/2)`, there is another option. If one provides after the name a value, the provided guess will be randomised with this factor. `0.02` would result into multiplying the guess with a random factor in [0.98,1.02]. The parameters which should be randomised are defined in the option `random_pars`. Note that if `random_guess` is not 0, the provided guess will be randomised in the selected parameters. If you want to use this option with the factor, set `random_guess` to 0.


Control file
============

+--------------+-------------------------------------------------------------------------------+
| Option       | Description                                                                   |
+==============+===============================================================================+
| cycle        | Number of cycles                                                              |
+--------------+-------------------------------------------------------------------------------+
| weights      | Weights for Stokes separated with `,`                                         |
+--------------+-------------------------------------------------------------------------------+
| nodes_temp   | Nodes for the temperature as list with `,`                                    |
+--------------+-------------------------------------------------------------------------------+
| nodes_magn   | Nodes for the magnetic field as list with `,`                                 |
+--------------+-------------------------------------------------------------------------------+
| nodes_vlos   | Nodes for the line-of-sight velocity as list with `,`                         |
+--------------+-------------------------------------------------------------------------------+
| nodes_gamma  | Nodes for the inclination of the magnetic field as list with `,`              |
+--------------+-------------------------------------------------------------------------------+
| nodes_phi    | Nodes for the azimuth of the magnetic field as list with `,`                  |
+--------------+-------------------------------------------------------------------------------+
| mu_cos       | Observer angle $\mu = \cos \theta$                                            |
+--------------+-------------------------------------------------------------------------------+
| abundace     | Abundance file                                                                |
+--------------+-------------------------------------------------------------------------------+
| gas_pressure | Gas Pressure at the top of the atmosphere                                     |
+--------------+-------------------------------------------------------------------------------+

Notes:
- For mode `2C`, the options nodes_* appear twice and each with a number (1 or 2) for each model

Radomisation Settings
=====================

+--------------+-------------------------------------------------------------------------------+
| Option       | Description                                                                   |
+==============+===============================================================================+
| random_guess | Number of random guesses (0 = use base model as guess)                        |
+--------------+-------------------------------------------------------------------------------+
| random_pars  | Parameters to be randomised separated as a list (e.g. `T,B,gamma,phi`)        |
+--------------+-------------------------------------------------------------------------------+
| lim_B        | Limits for the randomisation in B in G (e.g. `0,3000`)                        |
+--------------+-------------------------------------------------------------------------------+
| lim_vlos     |  Limits for the randomisation in v$_{los}$ in cm/s                            |
+--------------+-------------------------------------------------------------------------------+
| lim_gamma    |  Limits for the randomisation in the inclination in deg                       |
+--------------+-------------------------------------------------------------------------------+
| lim_phi      |  Limits for the randomisation in the azimuth in deg                           |
+--------------+-------------------------------------------------------------------------------+

Notes:
- For mode `2C`, the options lim_* appear twice and each with a number (1 or 2) for each model



---------
Mode `MC`
---------



General Stuff
=============

+--------------+----------------------------------------------------------------------------------+
| Option       | Description                                                                      |
+==============+==================================================================================+
| mode         | Mode (`MC`, `1C` or `2C`)                                                        |
+--------------+----------------------------------------------------------------------------------+
| path         | Path where everything is stored                                                  |
+--------------+----------------------------------------------------------------------------------+
| num          | Number of Models                                                                 |
+--------------+----------------------------------------------------------------------------------+
| model        | Model used as base model                                                         |
+--------------+----------------------------------------------------------------------------------+
| atoms        | Atoms used as a list where `;` = newline for the Grid File                       |
+--------------+----------------------------------------------------------------------------------+
| range_wave   | Range in the Grid (start in A, step in mA, number of wavelengths), `;` = newline |
+--------------+----------------------------------------------------------------------------------+

Data Stuff
==========

+--------------+-------------------------------------------------------------------------------+
| Option       | Description                                                                   |
+==============+===============================================================================+
| syn_out      | Prefix for the synthesis profiles and created models                          |
+--------------+-------------------------------------------------------------------------------+
| noise_out    | Prefix for the profiles with noise                                            |
+--------------+-------------------------------------------------------------------------------+
| inv_out      | Prefix for the output of the inversion results                                |
+--------------+-------------------------------------------------------------------------------+
| chi2         | Output of the chi2 values (npy)                                               |
+--------------+-------------------------------------------------------------------------------+


Creating Models and Synthesis
=============================

+--------------+-------------------------------------------------------------------------------+
| Option       | Description                                                                   |
+==============+===============================================================================+
| model_nodes  | Nodes for the models (options are 1, 2 or 3 nodes)                            |
+--------------+-------------------------------------------------------------------------------+
| model_pars   | Randomise these parameters for the created models                             |
+--------------+-------------------------------------------------------------------------------+
| noise_I      | Noise for Stokes I                                                            |
+--------------+-------------------------------------------------------------------------------+
| noise_Q      | Noise for Stokes Q                                                            |
+--------------+-------------------------------------------------------------------------------+
| noise_U      | Noise for Stokes U                                                            |
+--------------+-------------------------------------------------------------------------------+
| noise_V      | Noise for Stokes V                                                            |
+--------------+-------------------------------------------------------------------------------+
| create_B     | Limits for the first and last node in B (e.g. `0,4000;0,1000`)                |
+--------------+-------------------------------------------------------------------------------+
| create_vlos  | Limits for the first and last node in vlos (e.g. `-2e5,2e5;-2e5,2e5`)         |
+--------------+-------------------------------------------------------------------------------+
| create_gamma | Limits for the first and last node in gamma (e.g. `0,180;0,180`)              |
+--------------+-------------------------------------------------------------------------------+
| create_phi   | Limits for the first and last node in phi (e.g. `0,180;0,180`)                |
+--------------+-------------------------------------------------------------------------------+
| create_points| Limits are defined at these points (e.g. `1,-1,-4` for 3 nodes, '1,-4' for 2) |
+--------------+-------------------------------------------------------------------------------+

Note:
- The option `create_points` is not needed for `model_modes = 1`

Inversion configuration
=======================

+--------------+-------------------------------------------------------------------------------+
| Option       | Description                                                                   |
+==============+===============================================================================+
| line         | Lines File                                                                    |
+--------------+-------------------------------------------------------------------------------+
| guess        | Binary file used as initial guess, blank = base model                         |
+--------------+-------------------------------------------------------------------------------+
| cycle        | Number of cycles                                                              |
+--------------+-------------------------------------------------------------------------------+
| weights      | Weights for Stokes separated with `,`                                         |
+--------------+-------------------------------------------------------------------------------+
| nodes_temp   | Nodes for the temperature as list with `,`                                    |
+--------------+-------------------------------------------------------------------------------+
| nodes_magn   | Nodes for the magnetic field as list with `,`                                 |
+--------------+-------------------------------------------------------------------------------+
| nodes_vlos   | Nodes for the line-of-sight velocity as list with `,`                         |
+--------------+-------------------------------------------------------------------------------+
| nodes_gamma  | Nodes for the inclination of the magnetic field as list with `,`              |
+--------------+-------------------------------------------------------------------------------+
| nodes_phi    | Nodes for the azimuth of the magnetic field as list with `,`                  |
+--------------+-------------------------------------------------------------------------------+
| mu_cos       | Observer angle $\mu = \cos \theta$                                            |
+--------------+-------------------------------------------------------------------------------+
| abundace     | Abundance file                                                                |
+--------------+-------------------------------------------------------------------------------+
| gas_pressure | Gas Pressure at the top of the atmosphere                                     |
+--------------+-------------------------------------------------------------------------------+

Randomisation Settings
======================

+--------------+-------------------------------------------------------------------------------+
| Option       | Description                                                                   |
+==============+===============================================================================+
| random_guess | Number of random guesses (0 = use base model as guess)                        |
+--------------+-------------------------------------------------------------------------------+
| random_pars  | Parameters to be randomised separated as a list (e.g. `T,B,gamma,phi`)        |
+--------------+-------------------------------------------------------------------------------+
| lim_B        | Limits for the randomisation in B in G (e.g. `0,3000`)                        |
+--------------+-------------------------------------------------------------------------------+
| lim_vlos     |  Limits for the randomisation in v$_{los}$ in cm/s                            |
+--------------+-------------------------------------------------------------------------------+
| lim_gamma    |  Limits for the randomisation in the inclination in deg                       |
+--------------+-------------------------------------------------------------------------------+
| lim_phi      |  Limits for the randomisation in the azimuth in deg                           |
+--------------+-------------------------------------------------------------------------------+


