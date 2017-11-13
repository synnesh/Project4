import numpy as np
import matplotlib.pyplot as plt

import sys


info = np.loadtxt(sys.argv[1])

t=np.array(info[:,:1])
e_mean=np.array(info[:,1:2])
e_variance=np.array(info[:,2:3])
m_mean=np.array(info[:,3:4])
m_variance=np.array(info[:,4:5])
m_abs=np.array(info[:,5:6])
"""
accepted = np.array(info[:,6:7])
energy = np.array(info[:,7:8])

#mc=np.linspace(0,e_mean.size*1000, e_mean.size)
plt.figure()
data = [go.Histogram(energy=energy,histnorm='probability')]
py.iplot(data, filename='normalized histogram')
plt.title("Probability distribution for energy, variance = %0.4f" %e_variance[-1])
plt.xlabel("Energy states")
plt.ylabel("Probability")
plt.show()

#print mc.shape, e_mean.shape
"""
plt.figure(1)
plt.subplot(221)
plt.plot(t,e_mean)
plt.title("Mean energy as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Mean energy")

plt.subplot(222)
plt.plot(t,e_variance)
plt.title("Specific heat as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Specific heat")

plt.subplot(223)
plt.plot(t,m_mean)
plt.title("Magnetisation as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Magnetisation")

plt.subplot(224)
plt.plot(t,m_variance)
plt.title("Susceptibility as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")
plt.show()
"""
"""
plt.figure(2)
plt.plot(mc, e_mean)
plt.title("Mean energy as function of Monte Carlo cycles")
plt.grid()
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean energy")


plt.figure(3)
plt.plot(mc, m_mean)
plt.title("Mean magnetisation as function of Monte Carlo cycles")
plt.grid()
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean magnetisation")

plt.figure(3)
plt.plot(mc, m_mean)
plt.title("Mean magnetisation as function of Monte Carlo cycles")
plt.grid()
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean magnetisation")

plt.show()
