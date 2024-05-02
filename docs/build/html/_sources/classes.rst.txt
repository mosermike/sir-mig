Classes
=======

File model.py
-------------

Class model with all the physical parameter

.. autoclass:: model.Model
   :members: 

   .. automethod:: __init__
   .. automethod:: correct_phi
   .. automethod:: cut_to_map
   .. automethod:: get_attribute
   .. automethod:: interp
   .. automethod:: read
   .. automethod:: read_mod
   .. automethod:: read_results
   .. automethod:: set_dim
   .. automethod:: set_limit
   .. automethod:: write
   .. automethod:: write_model


File profile_stk.py
-------------------

Class profile_stk with the four Stokes Parameter

.. autoclass:: profile_stk.Profile
   :members:

   .. automethod:: __init__
   .. automethod:: __read_grid
   .. automethod:: __read_profile_sir
   .. automethod:: __read_profile_sir_mc
   .. automethod:: cut_to_map
   .. automethod:: cut_to_wave
   .. automethod:: read
   .. automethod:: read_results
   .. automethod:: read_results_MC
   .. automethod:: set_dim
   .. automethod:: write
   .. automethod:: write_profile
   .. automethod:: write_profile_mc
