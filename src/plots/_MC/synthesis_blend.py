"""
Plots the result of the SIR synthesis by blending two results into each other
"""

import numpy as np 
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../..")) 
import sir
import definitions as d
import matplotlib.pyplot as plt
from os.path import exists

print("NEEDS TO BE REVISED!")
sys.exit()
# Import library
# Import library
dirname = os.path.split(os.path.abspath(__file__))[0]
if exists(dirname + '/../../mml.mplstyle'):
	plt.style.use(dirname + '/../../mml.mplstyle')
elif "mml" in plt.style.available:
	plt.style.use('mml')



if sys.argv[1] == "-h" or (len(sys.argv) < 3):
     print("sir_synthesis_plot - Plots the result of a synthesis")
     print("Usage: python3 sir_synthesis_plot_blend [OPTION]")
     print()
     print("1: Results for model 1")
     print("2: Results for model 2")
     print("3: alpha (float or array)")
     print("4: Noise level as [Q, U, V] or one value for Q, U and V(optional, default no noise)")
     print("5: Additional save path (optional)")
     print("6: Additional text in filenames (optional)")
     print("7: Additional text in labels (optional)")
     if (len(sys.argv) < 3):
          print("[ERROR] Not enough arguments passed!")
     sys.exit()

# Define variables from input
Result1  = sys.argv[1] 			# Result for model 1
Result2  = sys.argv[2] 			# Result for model 1

# Alpha for defining the blending
if "[" in sys.argv[3]:
     alpha = sys.argv[3].replace(", ", ",").strip('][').split(',')
     alpha = [float(i) for i in alpha]
else:
     alpha  = [float(sys.argv[3])]

if len(sys.argv) > 4 and "[" in sys.argv[4]:
     noise = sys.argv[4].replace(", ", ",").strip('][').split(',')
     noise = [float(i) for i in noise]

elif len(sys.argv) > 4:
     noise = [float(sys.argv[4]), float(sys.argv[4]), float(sys.argv[4])]
else:
     noise = [0.0, 0.0, 0.0]

savepath = ''
if len(sys.argv) > 5:
     savepath = sys.argv[5]

add_text = ''
if len(sys.argv) > 6:
     add_text = sys.argv[6]

add_label = '_'
if len(sys.argv) > 7:
     add_label = sys.argv[7]

#############
# LOAD DATA #
#############
ll1, I1, Q1, U1, V1 = sir.read_profile(Result1)
ll2, I2, Q2, U2, V2 = sir.read_profile(Result2)

I = np.empty((len(alpha), I1.shape[0]))
Q = np.empty((len(alpha), Q1.shape[0]))
U = np.empty((len(alpha), U1.shape[0]))
V = np.empty((len(alpha), V1.shape[0]))

########################################################################
#                   PERFORM BLENDING                                   #
########################################################################
# Creating data which contain alpha percent of model 2 and (1-alpha) of model 1
for i in range(len(alpha)):
     I[i] = I1 * (1-alpha[i]) + I2 * alpha[i]
     Q[i] = Q1 * (1-alpha[i]) + Q2 * alpha[i]
     U[i] = U1 * (1-alpha[i]) + U2 * alpha[i]
     V[i] = V1 * (1-alpha[i]) + V2 * alpha[i]


########################################################################
#                   ADDING NOISE                                       #
########################################################################
if not noise[0] == noise[1] == noise[2] == 0.0:
     print("[STATUS] Noise added to Q, U and V")
     noise_Q = np.random.normal(scale = noise[0], size=len(Q[0]))
     noise_U = np.random.normal(scale = noise[1], size=len(U[0]))
     noise_V = np.random.normal(scale = noise[2], size=len(V[0]))

     for i in range(len(alpha)):
          Q[i] += noise_Q
          U[i] += noise_U
          V[i] += noise_V

     add_label += "\n with noise in Q,U,V"


fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12))
Markers = ["-", '--', 'dotted', 'dashdotdotted', 'densely dashed']

linestyle_str = [
     'solid',      # Same as (0, ()) or '-'
     'dotted',    # Same as (0, (1, 1)) or ':'
     'dashed',    # Same as '--'
     'dashdot',
     'solid',      # Same as (0, ()) or '-'
     'dotted',    # Same as (0, (1, 1)) or ':'
     'dashed',    # Same as '--'
     'dashdot',
     (0, (3,10,1,10))
     ] 

# Additional text in legend
#ax2.scatter([], [], color="w", alpha=0, label="Model 1 @ 2.5 kG")
#ax2.scatter([], [], color="w", alpha=0, label="Model 2 @ 8.0 kG")
ax2.scatter([], [], color="w", alpha=0, label=add_label)



ax1.plot(ll1, I1, label="Model 1 @ 2.5 kG", linewidth=3.0)
ax2.plot(ll1, Q1, label="Model 1 @ 2.5 kG", linewidth=3.0)
ax3.plot(ll1, U1, label="Model 1 @ 2.5 kG", linewidth=3.0)
ax4.plot(ll1, V1, label="Model 1 @ 2.5 kG", linewidth=3.0)

ax1.plot(ll2, I2, label="Model 2 @ 8.0 kG", linewidth=3.0)
ax2.plot(ll2, Q2, label="Model 2 @ 8.0 kG", linewidth=3.0)
ax3.plot(ll2, U2, label="Model 2 @ 8.0 kG", linewidth=3.0)
ax4.plot(ll2, V2, label="Model 2 @ 8.0 kG", linewidth=3.0)


for i in range(len(alpha)):
     # Delete first row as it is not important
     ax1.plot(ll1, I[i], label=r"$\alpha = $" + str(round(alpha[i]*100,1)) + r"$\%$",
              linestyle=linestyle_str[i])

     ax2.plot(ll1, Q[i], label=r"$\alpha = $" + str(round(alpha[i]*100,1)) + r"$\%$",
              linestyle=linestyle_str[i])

     ax3.plot(ll1, U[i], label=r"$\alpha = $" + str(round(alpha[i]*100,1)) + r"$\%$",
              linestyle=linestyle_str[i])

     ax4.plot(ll1, V[i], label=r"$\alpha = $" + str(round(alpha[i]*100,1)) + r"$\%$",
              linestyle=linestyle_str[i])

ax1.set_xlim(ll1[0], ll1[-1])
ax2.set_xlim(ll1[0], ll1[-1])
ax3.set_xlim(ll1[0], ll1[-1])
ax4.set_xlim(ll1[0], ll1[-1])

ax3.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')
ax4.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')

ax1.set_ylabel(r'$I/I_c$')
ax2.set_ylabel(r'$Q/I_c$')
ax3.set_ylabel(r'$U/I_c$')
ax4.set_ylabel(r'$V/I_c$')

ax2.legend(bbox_to_anchor=(1.01,0.95))

#plt.title("Blending of two synthesis")
# set the spacing between subplots
plt.tight_layout(pad=2.5)

plt.savefig(savepath + "sir_blended" + add_text)


