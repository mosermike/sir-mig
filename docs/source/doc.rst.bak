.. image:: pictures/sir-mig.png
   :alt: Icon of SIR-MIG
   :scale: 50%

Documentation
=============
 
Configuration
-------------

The code is always executed with a configuration file which is different for each mode. The options can be read here :doc:`here<config>`

Monte-Carlo Simulation
----------------------

The simulation is constructed the following way:

.. image:: pictures/mc.png
   :alt: Monte-Carlo-Simulation-Structure
   :scale: 50%

This can be done for as many numbers as needed. The code will do the following
     1. Create Models as described in the section //Constructing Models//.
     2. Synthesise the models to produce Stokes Profiles.
     3. Add Noise to the profiles.
     4. Invert the noisy profiles to get an Inversion Model.

Constructing Models
-------------------
Information about how the models are constructed can be found here:

- Creating Temperature Atmospheres: :doc:`temperature`

- Randomise Parameters: :doc:`Randomising Parameters<randomising>`

Adding Noise
------------

Noise is added by assuming a normal distribution with a given width $\sigma$ for each Stokes Parameter.

Inversion
---------

Contents
--------
.. toctree::
  :maxdepth: 0

  usage
  config
  randomising
  temperature
