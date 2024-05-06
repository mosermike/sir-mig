"""
Plots the changed parameters for one model
"""
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 
import sir
import definitions as d



# Import library
dirname = os.path.split(os.path.abspath(__file__))[0]
plt.rcParams["savefig.format"] = "pdf"
if d.plt_lib != "":
	plt.style.use(d.plt_lib)
else:
	if os.path.exists(dirname + '/mml.mplstyle'):
		plt.style.use(dirname + '/mml.mplstyle')
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
     print("sir_plot_one_model - Plot one model")
     print("Usage: python3 sir_plot_one_model [OPTION]")
     print()
     sir.option("[1. Pos]","Input model file")
     sir.option("-save","Savepath for plots (optional)")
     sir.option("-title","Title of the 4 figures plot")
     sir.option("-add","Additional text in filenames (optional)")
     sys.exit()


log_tau, T, Pe, vmicro, B, vlos,inc,azi,z, Pg, rho = sir.read_model(sys.argv[1])

# Additional labels
savepath = ''
if "-save" in sys.argv:
     savepath = sys.argv[sys.argv.index("-save")+1]

# Additional text in saved figures
add = ''
if "-add" in sys.argv:
     add = sys.argv[sys.argv.index("-add")+1]

# Title of the last plot
title = 'Example random model'
if "-title" in sys.argv:
     title = sys.argv[sys.argv.index("-title")+1]


###############################
#	Plot the 4 figures plot	#
###############################
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(16,12), sharex=True)

ax1.plot(log_tau, T)
ax2.plot(log_tau, B)
ax3.plot(log_tau, vlos/1e5)
ax4.plot(log_tau, inc)

################
#	Limits	#
################
ax1.set_xlim(log_tau[0], log_tau[-1])
ax2.set_xlim(log_tau[0], log_tau[-1])
ax3.set_xlim(log_tau[0], log_tau[-1])
ax4.set_xlim(log_tau[0], log_tau[-1])

################
#	Labels	#
################
ax4.set_xlabel(r"$\log \tau_{500}$")
ax3.set_xlabel(r"$\log \tau_{500}$")

titles   = [r"Temperature T", "","", r"Magnetic Field Strength B",
            r"Line-of-Sight Velocity $v_{\mathrm{los}}$", r"Inclination $\gamma$",
            r"Azimuth $\phi$"]

labels = ["T [K]", "", "", "B [G]", r"$v_{\mathrm{los}}$ $\left[\frac{\mathrm{km}}{\mathrm{s}}\right]$", r"$\gamma$ [deg]", r"$\phi$ [deg]"]

ax1.set_ylabel(labels[0])
ax2.set_ylabel(labels[3])
ax3.set_ylabel(labels[4])
ax4.set_ylabel(labels[5])


ax1.set_title(titles[0])
ax2.set_title(titles[3])
ax3.set_title(titles[4])
ax4.set_title(titles[5])

fig.suptitle(title, y=0.98, x=0.55)

plt.tight_layout()
plt.savefig(savepath + "one_model_pars" + add)

