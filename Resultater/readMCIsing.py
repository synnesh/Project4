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


