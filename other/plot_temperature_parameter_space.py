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

def rot(x,y,origin,deg):
	# Change the origin to the selected rotation point
	x_ = x - origin
	y_ = y - y[np.argmin(abs(x-origin))]

	# Normalise the functions
	x_max = np.max(np.abs(x_))
	y_max = np.max(np.abs(y_))
	x_ /= x_max
	y_ /= y_max

	# Rotate it by use of the rotation matrix
	rad = deg/180*np.pi
	x_rot =  x_*np.cos(rad)+y_*np.sin(rad)
	y_rot = -x_*np.sin(rad)+y_*np.cos(rad)

	# "Unormalise"
	x_rot *= x_max
	y_rot *= y_max

	# Interpolate in the selected temperature
	return y_rot+y[np.argmin(abs(x-origin))] #np.interp(x,np.flip(x_rot+origin),np.flip(y_rot+y[np.argmin(abs(x-origin))]))

def rot_ext(x1,y1,origin,deg1,deg2, func):
	T = np.zeros((2,len(x1)))
	x = np.copy(x1)
	y = np.copy(y1)
	T[0] = rot(x,y,origin,deg1)
	
	x = np.copy(x1)
	y = np.copy(y1)
	T[1] = rot(x,y,origin,deg2)
	

	for i in range(len(y1)):
			y1[i] = func([T[0][i],T[1][i]])

	return y1

deg_max = 3

# Compute the limits taking into account the rotation
T_maxs = []
T_max = rot_ext(log_taus,T_max,d.rot_point,-deg_max,deg_max, np.max)
for i in range(len(d.temp_f)):
	T_maxs.append(rot_ext(log_taus, cool_T +  d.temp_f[i][1]* HSRA_T,d.rot_point,-deg_max,deg_max,np.max))


T_min = rot_ext(log_taus,T_min,d.rot_point,-deg_max,deg_max, np.min)
T_mins = []
for i in range(len(d.temp_f)):
	T_mins.append(rot_ext(log_taus,cool_T +  d.temp_f[i][1]* HSRA_T,d.rot_point,-deg_max,deg_max,np.min))

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

def random_T(B):
	log_taus = np.copy(d.log_taus)
	HSRA_T = np.copy(d.upper_T)
	cool11_T = np.copy(d.lower_T)

	HSRA_T -= cool11_T  # Make it relative to cool11 so that with a fac of 1, I get the HSRA model

	# Factor for adding hsra depending on the magnetic field
	factor = np.random.uniform(d.lower_f,d.upper_f)
	# Apply restrictions for stronger magnetic fields
	for i in range(len(d.temp_B)):
		if B > d.temp_B[i]:
			factor = np.random.uniform(d.temp_f[i][0], d.temp_f[i][1])
			break
			
	# Little perturbation for cool11 model
	cool11_T = cool11_T * np.random.uniform(1-d.multiplicative_T, 1+d.multiplicative_T)
	# Add the values
	Ts = cool11_T + factor * HSRA_T
			
	# Rotate the temperature
	deg = np.random.uniform(-deg_max,deg_max)
	deg = deg_max
	Ts = rot(log_taus,Ts,d.rot_point,deg)
	
	return log_taus, Ts


fig, ax = plt.subplots()
#ax.set_prop_cycle(color=['0C5DA5', 'FF2C00', 'FF9500', '845B97', '474747', 'C20078','054907'])
#ax.plot(log_taus, d.upper_T, label="HSRA Model")
#ax.plot(log_taus, d.lower_T, label="Cool Umbra Model")
x1, y1 = random_T(0)
ax.plot(x1, y1, label="Random T")
x1, y1 = random_T(1000)
ax.plot(x1, y1, label="Random T")
x1, y1 = random_T(2500)
ax.plot(x1, y1, label="Random T")
x1, y1 = random_T(3500)
ax.plot(x1, y1, label="Random T")
x1, y1 = random_T(4000)
ax.plot(x1, y1, label="Random T")
x1, y1 = random_T(5000)
ax.plot(x1, y1, label="Random T")
ax.fill_between(log_taus, T_min,T_max, color="green",label="Parameter Space", alpha=0.5)
ax.set_xlim(log_taus[0],log_taus[-1])
ax.set_xlabel(r"$\log\tau_c$")
ax.set_ylabel(r"$T$ [K]")
ax.legend()
ax.set_title("Exemplary Random Temperature Profiles")
fig.savefig("temp_examples.pdf")


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

