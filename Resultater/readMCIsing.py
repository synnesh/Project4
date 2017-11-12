import numpy as np
import matplotlib.pyplot as plt

import sys

info = np.loadtxt(sys.argv[1])

t=np.transpose(info[:,:1])
e_mean=np.array(info[:,1:2])
e_variance=np.transpose(info[:,2:3])
m_mean=np.array(info[:,3:4])
m_variance=np.transpose(info[:,4:5])
m_abs=np.transpose(info[:,5:6])
accepted = np.array(info[:,6:7])
#print accepted
mc=np.linspace(0,e_mean.size*1000, e_mean.size)


#print mc.shape, e_mean.shape
plt.figure(1)
plt.plot(mc, e_mean)
plt.title("Mean energy as function of Monte Carlo cycles")
plt.grid()
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean energy")


plt.figure(2)
plt.plot(mc, m_mean)
plt.title("Mean magnetisation as function of Monte Carlo cycles")
plt.grid()
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean magnetisation")
plt.show()
