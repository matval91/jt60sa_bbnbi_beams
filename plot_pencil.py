import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

#=====================================================================================
# SET TEXT FONT AND SIZE
#=====================================================================================
plt.rc('font', family='serif', serif='Palatino')
#plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
#=====================================================================================


#VESSEL
vess_fname='2D_vess.txt'
vess=np.loadtxt(vess_fname, dtype=float, unpack=True)
R_vess, Z_vess = vess[0,:]*1e3, vess[1,:]*1e3
#---------------------------------------------------------------
fname='PNB_BBNBI.dat'
#fname='NNB_BBNBI.dat'
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

#figxy = plt.figure()
figxz = plt.figure()
#fig3d = plt.figure()
#axxy = figxy.add_subplot(111)
axxz = figxz.add_subplot(111)
#ax3d = fig3d.add_subplot(111, projection='3d')
##3D plot
#TOKAMAK
#shape of tokamak
phi = np.arange(0,2*np.pi,0.02*np.pi)
the = phi
x_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
y_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
for i,R in enumerate(R_vess):
    x_tok[:,i] = R*np.cos(phi)
    y_tok[:,i] = R*np.sin(phi)
z_tok = Z_vess

#pos_rect = [43,67,953,977]
pos_rect = [0, 1019]
pos_dir = 509

for i in range(len(data)):
    flag=0 #1->BB, 2->pencil, 3->narrow
    if fname == 'PNB_BBNBI.dat':
        if i<1020*1:
            flag=1
        #elif i<1020*13 and i>1020*12:
        #    flag=1
    else:
        if i%240==0:
            flag=1
    if flag==1:
        axxz.plot(pos[i,0], pos[i,2], 'k.')

    #if i in pos_rect: # Narrow beam
    #    flag=3

    if i == pos_dir:
        flag=2 #pencil model

    if flag==111:
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
        #axxy.plot([x1,x2],[y1,y2])
        #axxy.plot(x1,y1,'x')
        axxz.plot(x1,z1,'x')
        axxz.plot([x1,x2],[z1,z2])

    elif flag==2:
        [x1,y1,z1] = pos[pos_dir,:]
        print x1,y1,z1
        m=np.linalg.norm([x1,y1,z1])
        #m=20000
        #alpha is the one on xy plane, beta on xz plane
        ang_xy = math.atan2(y1,x1)+math.pi
        #print ang_xy, ang[i,0]
        alpha = ang_xy+ang[pos_dir,0]
        #print alpha
        #alpha = ang[i,0]
        beta  = ang[pos_dir,1]
        x2 = x1+math.cos(alpha)*math.cos(beta)*m
        y2 = y1+math.sin(alpha)*math.cos(beta)*m
        z2 = z1+math.sin(-beta)*m
        #axxy.plot([x1,x2],[y1,y2])
        #axxy.plot(x1,y1,'x')
        axxz.plot(x1,z1,'x')
        axxz.plot([x1,x2],[z1,z2], linewidth=2.5, color='k')

    elif flag==3:
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
        #axxy.plot([x1,x2],[y1,y2])
        #axxy.plot(x1,y1,'x')
        axxz.plot(x1,z1,'x')
        axxz.plot([x1,x2],[z1,z2], linewidth=2.5, color='k')

axxz.plot(R_vess, Z_vess,'m')
axxz.plot(-R_vess,Z_vess,'m')
#axxz.set_xlabel(r'R(mm)')
#axxz.set_ylabel(r'Z(mm)')
axxz.set_title(r'Narrow beam model')
axxz.set_xlim([-9500,-2000])
axxz.set_ylim([-7000,0])
axxz.tick_params(axis='x',which='both', bottom='off', top='off',labelbottom='off')
axxz.tick_params(axis='y',which='both', left='off', right='off',labelleft='off')

#axxy.plot(np.min(R_vess)*np.cos(the), np.min(R_vess)*np.sin(the),'m')
#axxy.plot(np.max(R_vess)*np.cos(the), np.max(R_vess)*np.sin(the),'m')
#ax3d.plot_surface(x_tok,y_tok,z_tok,color='k',alpha=0.2)

#axxy.scatter(intersections[:,0],intersections[:,1], color='y')
#axxz.scatter(intersections[:,0],intersections[:,2], color='y')

#axxz.axis('equal')
#axxy.axis('equal')
#ax3d.axis('equal')

plt.show()
