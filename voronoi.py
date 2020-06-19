import numpy as np
import pylab as pl

'''
points = []
n = 50
for i in range(n):
    points.append((np.random.rand(),np.random.rand()))

vp = Voronoi(points)
vp.process()
lines = vp.get_output()

px = []
py = []
for i in range(len(points)):
    px = np.hstack([px,points[i][0]])
    py = np.hstack([py,points[i][1]])


Z = []
ax = []
ay = []
bx = []
by = []
for i in range(len(lines)):
    ax = np.hstack([ax,lines[i][0]])
    ay = np.hstack([ay,lines[i][1]])
    bx = np.hstack([bx,lines[i][2]])
    by = np.hstack([by,lines[i][3]])
    z = np.sqrt((ax[-1]-bx[-1])**2+(ay[-1]-by[-1])**2)
    if(z<np.sqrt(2)):
        Z = np.hstack([Z,z])
'''

import h5py

f = h5py.File('points.h5','r')
px = f['x'][:]
py = f['y'][:]
ax = f['x1'][:]
ay = f['y1'][:]
bx = f['x2'][:]
by = f['y2'][:]

f.close()

'''
Z = np.sqrt((ax-bx)**2+(ay-by)**2)
ind = np.where(Z>np.sqrt(2))
Z = np.delete(Z,ind)
ax = np.delete(ax,ind)
ay = np.delete(ay,ind)
bx = np.delete(bx,ind)
by = np.delete(by,ind)
'''

f = h5py.File('cell.h5','r')
# test domain
nX = f['X'][:]
nY = f['Y'][:]
segLen = f['segL'][:]
segOrgX = f['segX'][:]
segOrgY = f['segY'][:]
segAng = f['segA'][:]


f.close()

print(ax.size)

F = pl.figure()
f = F.add_subplot(121)

f.plot(px,py,'.')
for i in range(ax.size):
    f.plot([ax,bx],[ay,by],'-',color='k')

#f.plot(X,Y,'*',color=(1,0,0))
#f.plot(Xeg,Yeg,'o',color=(0,1,0))
#f.plot(px[10],py[10],'*',color=(0,1,1))

#print(Xeg)

#f.axis([0,1,0,1])
#f.set_aspect(1)

#f = F.add_subplot(122)

f.plot(nX,nY,'o',color=(1,0,0))
f.plot(segOrgX,segOrgY,'*',color=(0,1,0))

for i in range(len(segOrgX)):
    xDst = segOrgX[i]+segLen[i]*np.cos(segAng[i])
    yDst = segOrgY[i]+segLen[i]*np.sin(segAng[i])
    f.plot([segOrgX[i],xDst],[segOrgY[i],yDst],'-',color=(1,1,0))

f.axis([0,1,0,1])
f.set_aspect(1)

print(nX)

pl.show()
