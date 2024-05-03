Functions
=========

This list provides information for each function in each file.

sir.py
------

Library for repeating functions such as reading the config, writing SIR files.

.. autofunction:: sir.list_to_string
.. autofunction:: sir.option
.. autofunction:: sir.read_chi2
.. autofunction:: sir.read_chi2s
.. autofunction:: sir.read_grid
.. autofunction:: sir.read_line
.. autofunction:: sir.read_config
.. autofunction:: sir.read_control
.. autofunction:: sir.read_model
.. autofunction:: sir.read_profile
.. autofunction:: sir.write_config_1c
.. autofunction:: sir.write_config_2c
.. autofunction:: sir.write_config_mc
.. autofunction:: sir.write_config
.. autofunction:: sir.write_control_1c
.. autofunction:: sir.write_control_2c
.. autofunction:: sir.write_control_mc
.. autofunction:: sir.write_grid
.. autofunction:: sir.write_grid_mc
.. autofunction:: sir.write_model
.. autofunction:: sir.write_profile


create_config.py
----------------

Functions to create the config files for the different modes.

.. autofunction:: create_config.config_MC
.. autofunction:: create_config.config_1C
.. autofunction:: create_config.config_2C


create_profiles
---------------
Function to create profiles with synthesis and adding noise
.. autofunction:: create_profiles.add_noise
.. autofunction:: create_profiles.synthesis


inversion.py
------------

Functions to perform the inversion

.. autofunction:: inversion.execute_inversion_1c
.. autofunction:: inversion.execute_inversion_2c
.. autofunction:: inversion.inversion_1c
.. autofunction:: inversion.inversion_2c
.. autofunction:: inversion.inversion_mc
.. autofunction:: inversion.scatter_data
.. autofunction:: inversion.scatter_data_mc


main.py
-------

Main file to start the program

.. autofunction:: main.main


model.py
--------

Function to read a binary model file with all the physical parameter

.. autofunction:: model.read_model

misc.py
-------

Miscellaneous functions

.. autofunction:: misc.create_task_folder_list
.. autofunction:: misc.initial
.. autofunction:: misc.x_y_add_zeros


obs.py
------

Functions related to the observation


.. autofunction:: obs.read_profile
.. autofunction:: obs.write_psf


profile_stk.py
--------------

Function to read a binary profile file with the four Stokes Parameter

.. autofunction:: profile_stk.read_profile


