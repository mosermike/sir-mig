"""
Plots the result of the SIR synthesis with different profiles
"""
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
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
dirname = os.path.split(os.path.abspath(__file__))[0]
plt.rcParams["savefig.format"] = "pdf"
if d.plt_lib != "":
	plt.style.use(d.plt_lib)
else:
	if exists(dirname + '/../mml.mplstyle'):
		plt.style.use(dirname + '/../mml.mplstyle')
		# if dvipng is not installed, dont use latex
		import shutil
		if shutil.which('dvipng') is None:
			plt.rcParams["text.usetex"] = "False"
			plt.rcParams["font.family"] = 'sans-serif'
			plt.rcParams["mathtext.fontset"] = 'dejavuserif'
	elif "mml" in plt.style.available:
		plt.style.use('mml')
		# if dvipng is not installed, dont use latex
		import shutil
		if shutil.which('dvipng') is None:
			plt.rcParams["text.usetex"] = "False"
			plt.rcParams["font.family"] = 'sans-serif'
			plt.rcParams["mathtext.fontset"] = 'dejavuserif'


if "-h" in sys.argv or (len(sys.argv) < 2):
     print("sir_synthesis_plot - Plots the result of a synthesis")
     print("Usage: python sir_synthesis_plot [OPTION]")
     print()
     print("1: Config file")
     print("2: Number(s) of the model(s)/profile(s) to be plotted")

     print("-label:  Labels (optional)")
     print("-save:   Additional save path (optional, default './')")
     print("-add:    Additional text in filenames (optional)")
     print("-title:  Title of the 4 figures plot")
     print("-title1: Title of the Stokes profiles")
     print("-gris:   Print gris information")
     print("-hinode: Print Hinode information")
     if (len(sys.argv) < 3):
          print("[ERROR] Not enough arguments passed!")
     sys.exit()

# Define variables from input
conf = sir.read_config(sys.argv[1])
path = conf["path"]
syn_out = os.path.join(path,conf["syn_out"])
syn_in = os.path.join(path,conf["syn_in"])

if "," in sys.argv[2]:
	Result = sys.argv[2].replace(", ", ",").strip('][').split(',')
	Result = [syn_out + i  + ".per" for i in Result]
else:
	Result  = [syn_out + sys.argv[2] + ".per"] 			# Result

# Synthesised model
if "," in sys.argv[2]:
	mod_syn = sys.argv[2].replace(", ", ",").strip('][').split(',')
	mod_syn = [syn_in + i + ".mod"for i in mod_syn]
else:
	mod_syn  = [syn_in + sys.argv[2] + ".mod"]



Labels = ["_" for i in Result]
if "-label" in sys.argv:
     ind = sys.argv.index("-save")+1
     if "," in sys.argv[ind]:
          Labels = sys.argv[ind].replace(", ", ",").strip('][').split(',')
     else:
          Labels  = [sys.argv[ind]] 			# Result





# Additional labels
savepath = ''
if "-save" in sys.argv:
     savepath = sys.argv[sys.argv.index("-save")+1]


# Additional text in saved figures
add_text = ''
if "-add" in sys.argv:
     add_text = sys.argv[sys.argv.index("-add")+1]

# Title of the last plot
title = ''
if "-title" in sys.argv:
     title = sys.argv[sys.argv.index("-title")+1]

title1 = ''
if "-title1" in sys.argv:
     title1 = sys.argv[sys.argv.index("-title1")+1]

datas = []
for i in range(len(Result)):
     datas.append(np.loadtxt(Result[i]))

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
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

for i in range(len(datas)):
     # Delete first row as it is not important
     data = np.delete(datas[i], 0, axis=1)
     if Labels[i] == '10kG':
          linewidth=3.0
     else:
          linewidth=2.0

     ax1.plot(data[:,0]/1000, data[:,1], label=Labels[i],
              linestyle=linestyle_str[i % len(linestyle_str)], linewidth=linewidth)

     ax2.plot(data[:,0]/1000, data[:,2], label=Labels[i],
              linestyle=linestyle_str[i % len(linestyle_str)], linewidth=linewidth)

     ax3.plot(data[:,0]/1000, data[:,3], label=Labels[i],
              linestyle=linestyle_str[i % len(linestyle_str)], linewidth=linewidth)

     ax4.plot(data[:,0]/1000, data[:,4], label=Labels[i],
              linestyle=linestyle_str[i % len(linestyle_str)], linewidth=linewidth)

ax1.set_xlim(data[0,0]/1000, data[-1,0]/1000)
ax2.set_xlim(data[0,0]/1000, data[-1,0]/1000)
ax3.set_xlim(data[0,0]/1000, data[-1,0]/1000)
ax4.set_xlim(data[0,0]/1000, data[-1,0]/1000)

#ax1.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')
#ax2.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')
if "-gris" in sys.argv:
     ax3.set_xlabel(r'$\Delta \lambda - 15648.514$ $[\AA]$', loc='center')
     ax4.set_xlabel(r'$\Delta \lambda - 15648.514$ $[\AA]$', loc='center')
elif "-hinode" in sys.argv:
     ax3.set_xlabel(r'$\Delta \lambda - 6301.5012$ $[\AA]$', loc='center')
     ax4.set_xlabel(r'$\Delta \lambda - 6301.5012$ $[\AA]$', loc='center')
else:
     ax3.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')
     ax4.set_xlabel(r'$\Delta \lambda$ $[\AA]$', loc='center')

ax1.set_ylabel(r'$I/I_c$')
ax2.set_ylabel(r'$Q/I_c$')
ax3.set_ylabel(r'$U/I_c$')
ax4.set_ylabel(r'$V/I_c$')

if "-gris" in sys.argv:
     fig.suptitle("Fe I 15645 Å, 15647 Å, 15648 Å + 15652 Å", y=0.98, x=0.55)
elif "-hinode" in sys.argv:
     fig.suptitle("Fe I 6301 Å + Fe I 6302 Å", y=0.98, x=0.55)
if title1 != '':
     fig.suptitle(title1, y=0.98, x=0.55)


if "-label" in sys.argv:
     ax2.legend(bbox_to_anchor=(1.01,0.95))


# set the spacing between subplots
plt.tight_layout(pad=2.5)

plt.savefig(path + "/" + savepath + "sir_synthesis" + add_text)


##
#    Plot the model
syn = []
for i in mod_syn:
     log_tau_syn, T_syn, Pe_syn, vmicro_syn, B_syn, vlos_syn, inc_syn, azi_syn, z_syn, Pg_syn, rho_syn = sir.read_model(i)
     syn.append([log_tau_syn, T_syn, Pe_syn, vmicro_syn, B_syn, vlos_syn, inc_syn, azi_syn, z_syn, Pg_syn, rho_syn])

labels = ["", "T [K]", r"$P_e$ $\left[\frac{\mathrm{dyn}}{\mathrm{cm}^2}\right]$", r"$v_{\mathrm{micro}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", "B [G]",
          r"$v_{\mathrm{los}}$ $\left[\frac{\mathrm{cm}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]"]
titles   = [r"Temperature T", r"Electron pressure $P_e$",r"Microturbulence velocity $v_{\mathrm{micro}}$", r"Magnetic Field Strength B", r"Line-of-Sight Velocity $v_{\mathrm{los}}$",
            r"Inclination $\gamma$"]

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, sharex=True)




for i in range(len(syn)):
     ax1.plot(syn[i][0], syn[i][1], label=Labels[i])
     ax2.plot(syn[i][0], syn[i][4], label=Labels[i])
     ax3.plot(syn[i][0], syn[i][5], label=Labels[i])
     ax4.plot(syn[i][0], syn[i][6], label=Labels[i])






# Set limits
ax1.set_xlim(syn[0][0][0], syn[0][0][-1])
ax2.set_xlim(syn[0][0][0], syn[0][0][-1])
ax3.set_xlim(syn[0][0][0], syn[0][0][-1])
ax4.set_xlim(syn[0][0][0], syn[0][0][-1])

# Set labels
ax3.set_xlabel(r"$\log \tau_{500}$")
ax4.set_xlabel(r"$\log \tau_{500}$")

ax1.set_ylabel(labels[1])
ax2.set_ylabel(labels[4])
ax3.set_ylabel(labels[5])
ax4.set_ylabel(labels[6])

ax1.set_title(titles[0])
ax2.set_title(titles[3])
ax3.set_title(titles[4])
ax4.set_title(titles[5])

if "-label" in sys.argv:
     ax2.legend(bbox_to_anchor=(1.01,0.95))

if title != '':
     fig.suptitle(title, y=0.98, x=0.55)

plt.tight_layout()

plt.savefig(path + "/" + savepath + "sir_synthesis_model" + add_text)


