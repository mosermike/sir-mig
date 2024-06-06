Usage
=====

.. _installation:

Installation
------------

To use SIR-MIG, first install the necessary libraries as defined in the text file requirements.py. The code might also work for other versions as mentioned there but it is not tested.

The necessary files with the right version can be installed by executing any of the two following lines in the directory of SIR-MIG:

.. code-block:: console

   (.venv) $ pip install pipreqs
   (.venv) $ conda install --yes --file requirements.txt -c conda-forge




Starting from real data
-----------------------
Here, how to start with fit files from GRIS data are explained. First, download the data you want to analyse.


Then, create a config file by executing
.. code-block:: console

   (.venv) $ python src/create_config.py

where several questions will be asked to create a config file which defines what the code will do. Copy or create all necessary files for the sir inversion (sir.x, abundance file, lines file, model file)

A description for the different options in the configuration files can be found :doc:`here <config>`

Afterwards, perform the main program by executing
.. code-block:: console

   (.venv) $ mpirun -np [num] python src/main.py [config file] -dir [path to l1 data]

Use the '-h' flag to print a help for the script 'main.py'.


Starting an MC simulation
-------------------------
Then, create a config file by executing
.. code-block:: console

   (.venv) $ python src/create_config.py

where several questions will be asked to create a config file which defines what the code will do.

Afterwards, one can start the simulation by performing
.. code-block:: console

   (.venv) $ mpirun -np [num] python src/main.py [config file]

Use the '-h' flag to print a help for the script 'main.py'.

