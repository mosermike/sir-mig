"""
Changes all possible values by a factor where log(tau) = 0 is taken as the reference point. Useful for changing the model. NO RANDOMIZATION IS HERE IMPLEMENTED.
"""

import numpy as np 
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../lib"))
import sir
from os.path import exists

if "-h" in sys.argv or (len(sys.argv) < 2):
	print("change_model - Changes components of a model to defined values")
	print("Usage: python change_model.py [OPTION]")
	print("The new model starts at the given value after the flag (e.g. -T 8000 => Model starts at 8000 K).")
	print()
	print("1: Input model file")
	print("2: Output model file")
	print("-T:	 New temperature in K (adding)")
	print("-Pe:	New electron pressure (multiplying)")
	print("-vmicro: New microturbulence velocity (multiplying)")
	print("-B:	 New magentic field strength in Gauss (multiplying)")
	print("-vlos:   New line of sight velocity in km/s (multiplying)")
	print("-gamma:    New inclination in deg (starting at this value, additive)")
	print("-phi:    New azimuth in deg  (starting at this value, additive)")
	print("-z:	 New real height in km (multiplying)")
	print("-Pg:	New gas pressure (multiplying)")
	print("-rho:    New density (multiplying)")
	sys.exit()

# Define variables from input
Input  = sys.argv[1] 		# Input model
Output = sys.argv[2]		# Output model

# Check if outputs exists
if exists(Output):
	inp = input("[WARN] Output " + Output + " exists. Do you want to continue? y/[n] ")
	while (inp != 'y'):
		
		if inp == 'n' or inp == '':
			print('Aborted')
			sys.exit()
		elif inp != 'y':
			inp = input('Input was not y or n. Try again. ')


##############################
# READ NEW VALUES FROM INPUT #
##############################
if "-T" in sys.argv:
	T_n	  = float(sys.argv[sys.argv.index("-T")+1])
	print("[NOTE] T is changed")

if "-Pe" in sys.argv:
	Pe_n	 = float(sys.argv[sys.argv.index("-Pe")+1])
	print("[NOTE] Pe is changed")

if "-vmicro" in sys.argv:
	v_micro_n = float(sys.argv[sys.argv.index("-vmicro")+1])
	print("[NOTE] vmicro is changed")

if "-B" in sys.argv:
	B_n	  = float(sys.argv[sys.argv.index("-B")+1])
	print("[NOTE] B is changed")

if "-vlos" in sys.argv:
	vlos_n    = float(sys.argv[sys.argv.index("-vlos")+1])*1e5
	print("[NOTE] vlos is changed")

if "-gamma" in sys.argv:
	inc_n	= float(sys.argv[sys.argv.index("-inc")+1])
	print("[NOTE] inc is changed")

if "-phi" in sys.argv:
	azi_n = float(sys.argv[sys.argv.index("-azi")+1])
	print("[NOTE] azi is changed")

if "-z" in sys.argv:
	z_n	  = float(sys.argv[sys.argv.index("-z")+1])
	print("[NOTE] z is changed")

if "-Pg" in sys.argv:
	Pg_n	 = float(sys.argv[sys.argv.index("-Pg")+1])
	print("[NOTE] Pg is changed")

if "-rho" in sys.argv:
	rho_n	= float(sys.argv[sys.argv.index("-rho")+1])
	print("[NOTE] rho is changed")

if ("-T" not in sys.argv and "-Pe" not in sys.argv and "-vmicro" not in sys.argv and 
	"-B" not in sys.argv and "-vlos" not in sys.argv and "-gamma" not in sys.argv and
	"-phi" not in sys.argv and "-z" not in sys.argv and "-Pg" not in sys.argv and
	"-rho" not in sys.argv):
	print("[WARN] No parameter is changed")
# Load data
File   = np.loadtxt(Input, skiprows=1)
Header = np.loadtxt(Input, max_rows=1)

# Create header
Header = "   " + str(Header[0]) + "  " + str(Header[1]) + "  " + str(Header[2])

# Create arrays with all the columns as rows

File_T = File.transpose()
log_tau = File_T[0]
T	  = File_T[1]
Pe	 = File_T[2]
v_micro = File_T[3]
B	  = File_T[4]
vlos    = File_T[5]
inc	= File_T[6]
azimuth = File_T[7]
if len(File_T) > 8:
	z	  = File_T[8]
	Pg	 = File_T[9]
	rho	= File_T[10]
else:
	z = None
	Pg = None
	rho = None

i0 = 0 #np.where(log_tau==0.0)[0][0] # Index where log tau = 0

####################################
#	    NEW TEMPERATURE		#
####################################
if "T_n" in vars():
	T = T - (T[i0]-T_n) # (T_n / T[i0]) * T
	if any(T < 2500):
		print("[Note] Perform hard cut so that T >= 2500.")
		T[T < 2500] = 2500

####################################
#	    NEW ELECTRON PRESSURE    #
####################################
if "Pe_n" in vars():
	if Pe[i0] == 0:
		Pe = np.ones(len(Pe))*Pe_n
	else:
		Pe = (Pe_n / Pe[i0]) * Pe

####################################
#	    NEW MICROTURBULENCE	 #
####################################
if "v_micro_n" in vars():
	if v_micro[i0] == 0:
		v_micro = np.ones(len(v_micro))*v_micro_n
	else:
		v_micro = (v_micro_n / v_micro[i0])*v_micro

####################################
#	    NEW MAGNETIC FIELD	  #
####################################
if "B_n" in vars():
	if B[i0] == 0:
		B = np.ones(len(B))*B_n
	else:
		B = (B_n / B[i0]) * B

####################################
#	    NEW VLOS			  #
####################################
if "vlos_n" in vars():
	if vlos[i0] == 0:
		vlos = np.ones(len(vlos))*vlos_n
	else:
		vlos = (vlos_n / vlos[i0]) * vlos

####################################
#	    NEW INCLINATION		#
####################################
if "inc_n" in vars():
		inc = inc + (inc_n-inc[0])

####################################
#	    NEW AZIMUTH		    #
####################################
if "azi_n" in vars():
		azimuth = azimuth + (azi_n-azimuth[0])

if (len(File_T) > 8):
	####################################
	#	    NEW HEIGHT			#
	####################################
	if "z_n" in vars():
		if z[i0] == 0:
			z = np.ones(len(z))*z_n
		else:
			z = z_n / z[i0] * z

	####################################
	#	    NEW GAS PRESSURE	    #
	####################################
	if "Pg_n" in vars():
		if Pg[i0] == 0:
			Pg = np.ones(len(Pg))*Pg_n
		else:
			Pg = Pg_n / Pg[i0] * Pg

	####################################
	#	    NEW DENSITY		    #
	####################################
	if "rho_n" in vars():
		if rho[i0] == 0:
			rho = np.ones(len(rho))*rho_n
		else:
			rho = rho_n / rho[i0] * rho


sir.write_model(Output, Header, log_tau, T, Pe,
					 v_micro, B, vlos, inc, azimuth, z, Pg, rho)
