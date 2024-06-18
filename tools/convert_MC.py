import numpy as np
import sys
sys.path(sys.path[0]+"/../src")
import model_atm as m
# Convert old version to new version of SIR-MIG-MC of the profiles
# 1st argument: first npy file
# 2nd argument: second npy file
# 3rd argument: destination


data1 = np.load(sys.argv[1])
data2 = np.load(sys.argv[2])

l1 = int(sys.argv[1].split("_")[-1].replace(".npy",""))
l2 = int(sys.argv[2].split("_")[-1].replace(".npy",""))

num1 = np.ones(data1.shape[2])*l1
num2 = np.ones(data2.shape[2])*l2


new = m.Model(data1.shape[0],1, data1.shape[2]+data2.shape[2])

new.indx = np.append(num1,num2)

new.wave[:,1,0:data1.shape[2]] = data1[:,0,:]
new.wave[:,1,data1.shape[2]:] = data2[:,0,:]

new.stki[:,2,0:data1.shape[2]] = data1[:,1,:]
new.stki[:,2,data1.shape[2]:] = data2[:,1,:]

new.stkq[:,3,0:data1.shape[2]] = data1[:,2,:]
new.stkq[:,3,data1.shape[2]:] = data2[:,2,:]

new.stku[:,4,0:data1.shape[2]] = data1[:,3,:]
new.stku[:,4,data1.shape[2]:] = data2[:,3,:]

new.stkv[:,5,0:data1.shape[2]] = data1[:,4,:]
new.stkv[:,5,data1.shape[2]:] = data2[:,4,:]

new.write(sys.argv[3])
