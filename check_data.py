import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('figure', facecolor='white')

#VESSEL
vess_fname='2D_vess.txt'
vess_fname='2D_div.txt'
vess=np.loadtxt(vess_fname, dtype=float, unpack=True)
R_vess, Z_vess = vess[0,:]*1e3, vess[1,:]*1e3
#---------------------------------------------------------------
fname='PNB_BBNBI.dat'
fname='NNB_BBNBI.dat'
with open(fname) as f:
    data=f.readlines()
pos = np.zeros((len(data),3),dtype=float)
ang = np.zeros((len(data),2),dtype=float)

if fname == 'PNB_BBNBI.dat':
    #INTERSECTIONS POINT
    intersections=np.array([[-2879.98234797,  2607.28513468,  -736.76402298],
       [-2879.98234797,  2607.28513468,   736.76402298],
       [-2067.0667934 ,  3287.65506232,  -735.689972  ],
       [-2067.0667934 ,  3287.65506232,   735.689972  ],
       [ 2067.99010851,  3289.87067538,  -737.48488931],
       [ 2067.99010851,  3289.87067538,   737.48488931],
       [ 2749.90416963,  4728.81113392,  -116.62958205],
       [ 2749.90416963,  4728.81113392,   116.62958205],
       [ 2750.08151154, -4728.91060227,  -116.6608528 ],
       [ 2750.08151154, -4728.91060227,   116.6608528 ],
       [-2067.99010851, -3289.87067538,  -737.48488931],
       [-2067.99010851, -3289.87067538,   737.48488931],
       [-2879.98234797,  2607.28513468,  -736.76402298],
       [-2879.98234797,  2607.28513468,   736.76402298],
       [-2067.0667934 ,  3287.65506232,  -735.689972  ],
       [-2067.0667934 ,  3287.65506232,   735.689972  ],
       [ 2067.99010851,  3289.87067538,  -737.48488931],
       [ 2067.99010851,  3289.87067538,   737.48488931],
       [ 2749.90416963,  4728.81113392,  -116.62958205],
       [ 2749.90416963,  4728.81113392,   116.62958205],
       [ 2750.08151154, -4728.91060227,  -116.6608528 ],
       [ 2750.08151154, -4728.91060227,   116.6608528 ],
       [-2067.99010851, -3289.87067538,  -737.48488931],
       [-2067.99010851, -3289.87067538,   737.48488931]])
     ###################################################
else:
    intersections=np.array([[-520.9488910742957,-2560.5421363981586, -799.3815],
                            [-520.9488910742957,-2560.5421363981586, -300.6185]])
for i in range(len(data)):
    pos[i,:]=data[i].split()[0:3]
    ang[i,:]=data[i].split()[3:5]
#print pos[0,:], ang[0,:]

figxy = plt.figure()
figxz = plt.figure()
fig3d = plt.figure()
axxy = figxy.add_subplot(111)
axxz = figxz.add_subplot(111)
ax3d = fig3d.add_subplot(111, projection='3d')
##3D plot
#TOKAMAK
#shape of tokamak
phi = np.arange(0,2.02*np.pi,0.02*np.pi)
the = phi
x_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
y_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
for i,R in enumerate(R_vess):
    x_tok[:,i] = R*np.cos(phi)
    y_tok[:,i] = R*np.sin(phi)
z_tok = Z_vess

for i in range(len(data)):
    flag=0
    if fname == 'PNB_BBNBI.dat':
        if i%1020==0:
            flag=1
            i+=530
        # if i<1020*1:
        #     flag=1
        # elif i<1020*13 and i>1020*12:
        #   flag=1
    else:
        if i%216==0:
            flag=1
            i+=103

    if flag==1:
        [x1,y1,z1] = pos[i,:]
        m=np.linalg.norm([x1,y1,z1])
        #m=20000
        #alpha is the one on xy plane, beta on xz plane

        ang_xy = math.atan2(y1,x1)+math.pi
        #print ang_xy, ang[i,0]

        alpha = ang_xy+ang[i,0]
        #print alpha
        #alpha = ang[i,0]
        beta  = ang[i,1]

        x2 = x1+math.cos(alpha)*math.cos(beta)*m
        y2 = y1+math.sin(alpha)*math.cos(beta)*m
        z2 = z1+math.sin(-beta)*m
        
        x1*=1e-3
        y1*=1e-3
        z1*=1e-3
        x2*=1e-3
        y2*=1e-3
        z2*=1e-3
        #plot of connecting lines between grid and focus point
        axxy.plot([x1,x2],[y1,y2], 'k--')
        #axxy.plot(x1,y1,'x')
        axxz.plot([x1,x2],[z1,z2],'k--')
        ax3d.plot([x1,x2],[y1,y2],zs=[z1,z2], c='k')

#plot of vessell
axxz.plot(R_vess*1e-3, Z_vess*1e-3,'k')
axxz.plot(-R_vess*1e-3,Z_vess*1e-3,'k')
axxy.plot(np.min(R_vess)*np.cos(the)*1e-3, np.min(R_vess)*np.sin(the)*1e-3,'k')
axxy.plot(np.max(R_vess)*np.cos(the)*1e-3, np.max(R_vess)*np.sin(the)*1e-3,'k')
#ax3d.plot_surface(x_tok*1e-3,y_tok*1e-3,z_tok*1e-3,color='k',alpha=0.2)

#Plotting beamlets points
axxy.scatter(pos[:,0]*1e-3, pos[:,1]*1e-3, 20, 'k')
axxz.scatter(pos[:,0]*1e-3, pos[:,2]*1e-3, 10, 'k')
axxy.scatter(intersections[:,0]*1e-3,intersections[:,1]*1e-3, color='k')
axxz.scatter(intersections[:,0]*1e-3,intersections[:,2]*1e-3, color='k')
ax3d.scatter(pos[:,0]*1e-3, pos[:,1]*1e-3, zs=pos[:,2]*1e-3, s=20, c='k')
axxz.axis('equal')
axxy.axis('equal')
ax3d.set_aspect('equal')
#plt.setp(axxy.get_yticklabels()[0], visible=False) 
#plt.setp(axxz.get_yticklabels()[0], visible=False) 
#plt.setp(ax3d.get_yticklabels()[0], visible=False) 

axxy.set_xlabel(r'X [m]'); axxy.set_ylabel(r'Y [m]')
axxz.set_xlabel(r'R [m]'); axxz.set_ylabel(r'Z [m]')
ax3d.set_xlabel(r'X [m]'); ax3d.set_ylabel(r'Y [m]'); ax3d.set_zlabel(r'Z [m]')

axxy.grid('on')
axxz.grid('on')
ax3d.grid('on')
axxz.set_xlim([-10,10]); axxy.set_xlim([-10,10]); axxy.set_ylim([-10,10])

figxy.tight_layout()
figxz.tight_layout()
fig3d.tight_layout()
plt.show()
