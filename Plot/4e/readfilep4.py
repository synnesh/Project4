import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
from scipy.optimize import curve_fit


info1 = np.loadtxt(sys.argv[1])

t1=np.array(info1[:,:1])
e1_mean=np.array(info1[:,1:2])
e1_variance=np.array(info1[:,2:3])
m1_mean=np.array(info1[:,3:4])
m1_variance=np.array(info1[:,4:5])
m1_abs=np.array(info1[:,5:6])

info2 = np.loadtxt(sys.argv[2])

t2=np.array(info2[:,:1])
e2_mean=np.array(info2[:,1:2])
e2_variance=np.array(info2[:,2:3])
m2_mean=np.array(info2[:,3:4])
m2_variance=np.array(info2[:,4:5])
m2_abs=np.array(info2[:,5:6])

info3 = np.loadtxt(sys.argv[3])

t3=np.array(info3[:,:1])
e3_mean=np.array(info3[:,1:2])
e3_variance=np.array(info3[:,2:3])
m3_mean=np.array(info3[:,3:4])
m3_variance=np.array(info3[:,4:5])
m3_abs=np.array(info3[:,5:6])

info4 = np.loadtxt(sys.argv[4])

t4=np.array(info4[:,:1])
e4_mean=np.array(info4[:,1:2])
e4_variance=np.array(info4[:,2:3])
m4_mean=np.array(info4[:,3:4])
m4_variance=np.array(info4[:,4:5])
m4_abs=np.array(info4[:,5:6])
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

plt.plot(t1,e1_mean, t2, e2_mean, t3, e3_mean, t4, e4_mean)
plt.title("Mean energy as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Mean energy")
plt.legend(["L=40", "L=60", "L=80","L=100"], loc=2)
#plt.show()


plt.plot(t1,e1_variance, t2,e2_variance, t3,e3_variance, t4, e4_variance)
plt.title("Specific heat as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Specific heat")
plt.legend(["L=40", "L=60", "L=80","L=100"], loc=2)
#plt.show()

plt.plot(t1,m1_abs, t2,m2_abs, t3,m3_abs, t4,m4_abs)
plt.title("Magnetisation as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Magnetisation")
plt.legend(["L=40", "L=60", "L=80","L=100"], loc=2)
#plt.show()

plt.plot(t1,m1_variance, t2,m2_variance, t3,m3_variance, t4,m4_variance)
plt.title("Susceptibility as function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")
plt.legend(["L=40", "L=60", "L=80","L=100"], loc=2)


#plt.show()

susc100=interp1d(np.transpose(t4) [0],np.transpose(m4_variance)[0], kind='cubic')
susc80=interp1d(np.transpose(t3) [0],np.transpose(m3_variance)[0], kind='cubic')

specheat100=interp1d(np.transpose(t4) [0],np.transpose(e4_variance)[0], kind='cubic')
specheat80=interp1d(np.transpose(t3) [0],np.transpose(e3_variance)[0], kind='cubic')


tnew=np.linspace(2.15, 2.35, 10000)

fsusc100=susc100(tnew)
fsusc80=susc80(tnew)

fspecheat100=specheat100(tnew)
fspecheat80=specheat80(tnew)

su100index=0
su80index=0
sh100index=0
sh80index=0

su100=0
su80=0
sh100=0
sh80=0

for i in range(len(tnew)):	
	if(fspecheat100[i]>sh100):
		sh100=fspecheat100[i]
		sh100index=i
	if(fspecheat80[i]>sh80):
		sh80=fspecheat80[i]
		sh80index=i
	if(fsusc100[i]>su100):
		su100=fsusc100[i]
		su100index=i
	if(fsusc100[i]>su100):
		su100=fsusc100[i]
		su100index=i

print "susc100 temperature: %f" %tnew[su100index]
print "susc80 temperature: %f" %tnew[su80index]
print "sh100 temperature: %f" %tnew[sh100index]
print "sh80 temperature: %f" %tnew[sh80index]


