import numpy as np
import sys
sys.path.append(sys.path[0]+"/../")
import model_atm as m
# Convert old version to new version of SIR-MIG-MC of the models
# 1st argument: first file
# 2nd argument: destination


data = np.load(sys.argv[1])

mod = m.model_atm(data.shape[0],data.shape[1], data.shape[3])

mod.tau = data[0,0,0,:]
mod.T = data[:,:,1]
mod.Pe = data[:, :, 2]
mod.vmicro = data[:, :, 3]
mod.B = data[:, :, 4]
mod.vlos = data[:, :, 5]
mod.gamma = data[:, :, 6]
mod.phi = data[:, :, 7]
mod.z = data[:, :, 8]
mod.Pg = data[:, :, 9]
mod.rho = data[:, :, 10]
mod.fill += 1
mod.vmacro += 0.1

mod.write(sys.argv[2])
