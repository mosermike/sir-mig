import sys
sys.path.append("../src/")
sys.path.append("../plots/")
import definitions as d
import matplotlib.pyplot as plt
import sir
import numpy as np

plt.ion()
# Import library
sir.mpl_library()


d.multiplicative_T


log_taus = d.log_taus
HSRA_T = np.copy(d.upper_T)
cool_T = np.copy(d.lower_T)
		
HSRA_T -= cool_T # Make it relative to cool11 so that with a fac of 1, I get the HSRA model
		
# Little perturbation for cool model
T_min = cool_T * (1-d.multiplicative_T)
		
T_max = cool_T + d.upper_f * HSRA_T

# Add (only in creating models) additional perturbation in a resulting rotation around log tau -1
#factor = np.random.uniform(0.9, 1.1)
def rot_max(T,f1,f2):
	T_maxs = np.zeros(shape=(2,len(T)))

	Ts = np.copy(T)
	Ts[Ts > d.rot_point] = Ts[Ts > d.rot_point] * f1
	Ts[Ts <= d.rot_point] = Ts[Ts <= d.rot_point] / f1
	T_maxs[0] = Ts
	Ts = np.copy(T)
	Ts[Ts > d.rot_point] = Ts[Ts > d.rot_point] * f2
	Ts[Ts <= d.rot_point] = Ts[Ts <= d.rot_point] / f2
	T_maxs[1] = Ts

	for i in range(len(T)):
		T[i] = np.max([T_maxs[0][i],T_maxs[1][i]])
	return T

T_max = rot_max(T_max,1.05,0.95)

T_maxs = []
for i in range(len(d.temp_f)):
	T_maxs.append(rot_max(cool_T +  d.temp_f[i][1]* HSRA_T,1.05,0.95))

T_mins = np.zeros(shape=(2,len(T_min)))
Ts = np.copy(T_min)
Ts[Ts > d.rot_point] = Ts[Ts > d.rot_point] * 1.05
Ts[Ts <= d.rot_point] = Ts[Ts <= d.rot_point] / 1.05
T_mins[0] = Ts
Ts = np.copy(T_min)
Ts[Ts > d.rot_point] = Ts[Ts > d.rot_point] * 0.95
Ts[Ts <= d.rot_point] = Ts[Ts <= d.rot_point] / 0.95
T_mins[1] = Ts

for i in range(len(T_min)):
	T_min[i] = np.min([T_mins[0][i],T_mins[1][i]])


fig, ax = plt.subplots()
ax.set_prop_cycle(color=['0C5DA5', 'FF2C00', 'FF9500', '845B97', '474747', 'C20078','054907'])
ax.plot(log_taus, d.upper_T, label="HSRA Model")
ax.plot(log_taus, d.lower_T, label="Cool Umbra Model")
ax.fill_between(log_taus, T_min,T_max, color="green",label="Parameter Space", alpha=0.5)
ax.set_xlim(log_taus[0],log_taus[-1])
ax.set_xlabel(r"$\log\tau_c$")
ax.set_ylabel(r"$T$ [K]")
ax.legend()
ax.set_title("Creation of Random Temperature Profiles")
fig.savefig("temp_pars_space.pdf")

fig, ax = plt.subplots()
ax.set_prop_cycle(color=['0C5DA5', 'FF2C00', 'FF9500', '845B97', '474747', 'C20078','054907'])
ax.plot(log_taus, d.upper_T, label="HSRA Model")
ax.plot(log_taus, d.lower_T, label="Cool Umbra Model")
ax.fill_between(log_taus, T_min,T_max, color="green",label="Parameter Space", alpha=0.5)
for i in range(len(d.temp_f)):
	ax.plot(log_taus, T_maxs[i], linestyle="dotted", label=r"Upper Limit for $B$ $>$ " + str(d.temp_B[i]) + r" G", alpha=0.6)
ax.set_xlim(log_taus[0],log_taus[-1])
ax.set_xlabel(r"$\log\tau_c$")
ax.set_ylabel(r"$T$ [K]")
ax.legend()
ax.set_title("Creation of Random Temperature Profiles")
fig.savefig("temp_pars_space_upper_lim.pdf")



T_min = cool_T * (1-d.multiplicative_T)
		
T_max = cool_T + d.upper_f * HSRA_T

fig, ax = plt.subplots()
ax.set_prop_cycle(color=['0C5DA5', 'FF2C00', 'FF9500', '845B97', '474747', 'C20078','054907'])
ax.plot(log_taus, d.upper_T, label="HSRA Model")
ax.plot(log_taus, d.lower_T, label="Cool Umbra Model")
ax.fill_between(log_taus, T_min,T_max, color="green",label="Parameter Space", alpha=0.5)
ax.set_xlim(log_taus[0],log_taus[-1])
ax.set_xlabel(r"$\log\tau_c$")
ax.set_ylabel(r"$T$ [K]")
ax.legend()
ax.set_title("Creation of Random Temperature Profiles for Inversion")
fig.savefig("temp_pars_space_inversion.pdf")

fig, ax = plt.subplots()
ax.set_prop_cycle(color=['0C5DA5', 'FF2C00', 'FF9500', '845B97', '474747', 'C20078','054907'])
ax.plot(log_taus, d.upper_T, label="HSRA Model")
ax.plot(log_taus, d.lower_T, label="Cool Umbra Model")
ax.fill_between(log_taus, T_min,T_max, color="green",label="Parameter Space", alpha=0.5)
for i in range(len(d.temp_f)):
	ax.plot(log_taus, cool_T + d.temp_f[i][1] * HSRA_T, linestyle="dotted", label=r"Upper Limit for $B$ $>$ " + str(d.temp_B[i]) + r" G", alpha=0.6)
ax.set_xlim(log_taus[0],log_taus[-1])
ax.set_xlabel(r"$\log\tau_c$")
ax.set_ylabel(r"$T$ [K]")
ax.legend()
ax.set_title("Creation of Random Temperature Profiles for Inversion")
fig.savefig("temp_pars_space_upper_lim_inversion.pdf")

