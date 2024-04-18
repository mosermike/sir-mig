import numpy as np
import sys

# Convert old version to new version of SIR-MIG-MC
# 1st argument: first npy file
# 2nd argument: second npy file
# 3rd argument: destination


data1 = np.load(sys.argv[1])
data2 = np.load(sys.argv[2])

l1 = int(sys.argv[1].split("_")[-1].replace(".npy",""))
l2 = int(sys.argv[2].split("_")[-1].replace(".npy",""))

num1 = np.ones(data1.shape[2])*l1
num2 = np.ones(data2.shape[2])*l2


new = np.zeros(shape=(data1.shape[0],6,data1.shape[2]+data2.shape[2]))

new[:,0,:] = np.append(num1,num2)

new[:,1,0:data1.shape[2]] = data1[:,0,:]
new[:,1,data1.shape[2]:] = data2[:,0,:]

new[:,2,0:data1.shape[2]] = data1[:,1,:]
new[:,2,data1.shape[2]:] = data2[:,1,:]

new[:,3,0:data1.shape[2]] = data1[:,2,:]
new[:,3,data1.shape[2]:] = data2[:,2,:]

new[:,4,0:data1.shape[2]] = data1[:,3,:]
new[:,4,data1.shape[2]:] = data2[:,3,:]

new[:,5,0:data1.shape[2]] = data1[:,4,:]
new[:,5,data1.shape[2]:] = data2[:,4,:]

np.save(sys.argv[3], new)
