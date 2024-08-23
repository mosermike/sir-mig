.. image:: pictures/sir-mig.png
   :alt: Icon of SIR-MIG
   :scale: 50%


SIR-MIG Documentation
=====================

**Welcome to SIR MIG's documentation!**

.. toctree::
   :maxdepth: 1
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Description
===========

**SIR-MIG** is an algorithm based on the inversion code SIR. It parallelises SIR and implements functionalities to use multiple random guesses.
There are three different modes selectable:

- 1C : 1 Component Inversion
- 2C : 2 Component Inversion
- MC : Monte-Carlo Simulation
- SY : Synthesis

It is also possible to use it without random guesses.


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

.. note::

   This project is under active development.

Contents
--------
.. toctree::
  :maxdepth: 1

  usage
  config
  randomising
  temperature
  api
