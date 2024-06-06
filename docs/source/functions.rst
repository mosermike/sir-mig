Functions
=========

The code has several files with different functions in them. They all serve different purposes


This list provides information for each function in each file.

sir.py
------

+-----------------------------+---------------------------------------------------------------------------------------------+
| Function                    | Description                                                                                 |
+=============================+=============================================================================================+
| sir.create_task_folder_list | Creates a list which folders should be created and executed. This is done so that the       |
|                             | inversion itself can be executed linearly to make use of all cores.                         |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.initial                 | Initial print outs and preparation                                                          |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.list_to_string          | Convert a list to a string                                                                  |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.option                  | Prints out an option in the help page                                                       |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_chi2               | Reads the last chi value in a inv.chi file                                                  |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_chi2s              | Reads all the chi2 from the inversion                                                       |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_config             | Reads a config file for the inversion                                                       |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_control            | Reads a control file in the scheme SIR expects it.                                          |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_grid               | Reads the grid file                                                                         |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_line               | Reads the line file                                                                         |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_model              | Reads a model file and returns all parameters                                               |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.read_profile            | Reads the first LINE data from a profile computed by SIR                                    |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.write_config            | Writes a config file with the information provided as a dictionary                          |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.write_control           | Writes a control file in the scheme SIR expects it.                                         |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.write_gauss_psf         | Writes the spectral point spread function with the given sigma. This function is used when  |
|                             | the field psf in the config contains ‘gauss 1.0’ with sigma = 1.0 mA in this example.       |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.write_grid              | Writes the Grid file with data from the config file                                         |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.write_model             | Writes a model with the given data in a specific format. Note that negative values have one |
|                             | white space less                                                                            |
+-----------------------------+---------------------------------------------------------------------------------------------+
| sir.x_y_add_zeros(x, y)     | Adds zeros so that the returning strings have 4 letters                                     |
+-----------------------------+---------------------------------------------------------------------------------------------+

A full description of each function can be found :doc:`here<sir>`.

create_config.py
----------------

.. automodule:: create_config
  :members:
  :no-index:

model_atm.py
------------

Function to read a binary model file with all the physical parameter

.. autofunction:: model_atm.read_model
  :no-index:

profile_stk.py
--------------

Function to read a binary profile file with the four Stokes Parameter

.. autofunction:: profile_stk.read_profile
  :no-index:

Contents
--------
.. toctree::
  :maxdepth: 1

  sir
