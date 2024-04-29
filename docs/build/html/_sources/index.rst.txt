.. Temp documentation master file, created by
   sphinx-quickstart on Mon Apr 29 11:33:54 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SIR MIG's documentation!
===================================

.. toctree::
   :maxdepth: 2
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

It is also possible to use it without random guesses.

Check out the :doc:`usage` section for further information.

Check out the :doc:`functions` section for a list of functions.

Check out the :doc:`classes` section for a list of classes which are used to read and write the Stokes profiles and the models.


.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   functions
   classes
   api
   
