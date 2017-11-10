import numpy as np
import matplotlib.pyplot as plt7

import sys

info = np.loadtxt(sys.argv[1])

t=np.transpose(info[:,:1])
e_mean=np.transpose(info[:,1:2])
e_variance=np.transpose(info[:,2:3]
m_mean=np.transpose(info[:,3:4])
m_variance=np.transpose(info[:,4:5])
m_abs=np.transpose(info[:,5:6])

mc_cycles=[]
for i in range(len(e_mean)):
	mc_cycles.append(i*1000)
	
mc=np.array(mc_cycles)

plt.plot(mc, e_mean)
plt.show()
