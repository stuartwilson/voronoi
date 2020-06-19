import numpy as np
import pylab as pl
import h5py

F = h5py.File('ising.h5','r')
T = F['T'][:]
M = F['M'][:]
E = F['E'][:]
V = F['V'][:]
F.close()

#x = np.log(T)
x = T

F = pl.figure()
f = F.add_subplot(311)
f.plot(x,M,color=(0,0,0))
f.set_xlim(1,np.max(x))
f.set_xscale('log')
f.set_xlabel('Temperature')
f.set_ylabel('Magnetization')

f = F.add_subplot(312)
f.plot(x,E,color=(0,0,0))
for i in range(len(V)):
    f.plot([x[i],x[i]],[E[i]-V[i],E[i]+V[i]],color=(0,0,0))
f.set_xlim(1,np.max(x))
f.set_xscale('log')
f.set_xlabel('Temperature')
f.set_ylabel('Energy')

f = F.add_subplot(313)
f.plot(x,V,color=(0,0,0))
f.set_xlim(1,np.max(x))
f.set_xscale('log')
f.set_xlabel('Temperature')
f.set_ylabel('stdv. Energy')


pl.show()
