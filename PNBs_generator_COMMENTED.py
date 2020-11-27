# Script to generate the beamlet data for the JT60-SA positive NBs
# with data obtained from the references in the JT60-SA document
# M. Vallar - 05/2017
#-----------------------------------------------------------------

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

def rotate_point(x,y,theta):
    """
    Calculate the rotation of a point given the angle theta
    """
    xnew=x*math.cos(theta)-y*math.sin(theta)
    ynew=x*math.sin(theta)+y*math.cos(theta)
    return xnew, ynew

def find_basis(vector):
    """
    Find the basis for the plane defined by vector. One vector is in the form (a,b,0), w a*a+b*b=1
    """
    b=vector[0]/math.sqrt(vector[0]*vector[0]+vector[1]*vector[1])
    a=-b*vector[1]/vector[0]
    basis_1=[a,b,0]
    basis_2=np.cross(basis_1,vector)
    basis_2=basis_2/np.linalg.norm(basis_2)
    return basis_1, basis_2

def line_intersection(line1, line2):
    """
    http://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python
    """
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
    x,y = [0,0]
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        print "Lines do not intersect"
        return


    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


plot_flag=1

# DESIGN
#POSITIVE GRID: 12x27 cm2
tot_size = [float(268.8), float(117.3)] # total size of the grid in mm [row, col]
ap_diam = 18 # diameter of a beamlet aperture in mm
# DATA RELATIVE TO THE SEGMENTS
ap_grid = [43, 24] # number of apertures in segment [row, col]
# THE FIRST AND LAST LINE HAVE 4 HOLES LESS         (SO 20 BEAMLETS)
# THE SECOND AND BEFORE-LAST LINE HAVE 2 HOLES LESS (SO 22 BEAMLETS)
h_ap_dist = 5.1
v_ap_dist = 6.4
#-------------------------------------------

#VESSEL
vess_fname='2D_vess.txt'
vess=np.loadtxt(vess_fname, dtype=float, unpack=True)
R_vess, Z_vess = vess[0,:]*1e3, vess[1,:]*1e3
#---------------------------------------------------------------
#COORDINATES OF TF COIL1:
ymin_TFCOIL = -0.3*1000
ymax_TFCOIL = -ymin_TFCOIL
xmin_TFCOIL = 1.6*1000
xmax_TFCOIL = 4.7*1000
w_TF=3*1000
h_TF=1*1000
#ROTATION FOR CORRECT TOROIDAL 0 OF TF COIL
ang_rotation = -11*math.pi/180.  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ANGLE OF ROTATION (+: CCW, -:CW)
x1_TF,y1_TF = rotate_point(xmin_TFCOIL, ymin_TFCOIL, ang_rotation)
x2_TF,y2_TF = rotate_point(xmin_TFCOIL, ymax_TFCOIL, ang_rotation)
x3_TF,y3_TF = rotate_point(xmax_TFCOIL, ymax_TFCOIL, ang_rotation)
x4_TF,y4_TF = rotate_point(xmax_TFCOIL, ymin_TFCOIL, ang_rotation)
#---------------------------------------------------------------

# DATA RELATIVE TO THE INJECTORS FROM EXCEL FILE
data_fname = 'AB_fromexcel.dat'
names_label=['1A','2A','3A','4A','5A','6A','7A','8A','9A','10A',\
             '13A','14A','1B','2B','3B','4B','5B','6B','7B','8B',\
             '9B','10B','13B','14B']
data_f = open(data_fname)
data = np.loadtxt(data_f, dtype=float, skiprows=2, unpack=True)
[xs_a,ys_a,zs_a] = data[0:3,:]
[theta_a, phi_a] = data[3:5,:]*math.pi/180.
data_f.close()
#---------------------------------------------------------------


#ROTATION FOR CORRECT TOROIDAL 0 OF SOURCES
r_a = np.sqrt(xs_a*xs_a+ys_a*ys_a)
original_angle = np.arctan2(ys_a,xs_a)
xs_a = r_a*np.cos(original_angle+ang_rotation)
ys_a = r_a*np.sin(original_angle+ang_rotation)
#---------------------------------------------------------------
#ROTATION FOR CORRECT TOROIDAL 0 OF AIMING VECTORS
phi_a = phi_a-ang_rotation
#---------------------------------------------------------------

n_beams = xs_a.shape[0]
print "Number of beam sources (A/B): ", n_beams

# Initialisation of data
v_ort = np.zeros((n_beams, 3), dtype=float)
centers = np.zeros((n_beams, ap_grid[1], ap_grid[0], 3), dtype=float) #Centers x,y,z position
centers_new = np.copy(centers)
intersections = np.zeros((n_beams,3), dtype=float)
v_ort_bmlts = np.zeros((n_beams, ap_grid[1], ap_grid[0], 3), dtype=float)
# XYZ wrt to machine center
x_int, y_int, z_int = 0,0,0

# VECTOR NORMAL TO THE GRID JUST CREATED:
orth_vect = [1,0,0]
#VECTOR TO BE NORMAL TO THE GRID
aim_vect = np.zeros((n_beams,3), dtype=float) #These is the set of all the aiming vectors
for i in range(n_beams):
    #alpha is the one on xy plane, beta on xz plane
    alpha = -phi_a[i]
    beta  = theta_a[i]
    x = math.cos(alpha)*math.cos(beta)
    y = math.sin(alpha)*math.cos(beta)
    z = math.sin(beta)
    aim_vect[i,:] = [x,y,z]


# vector from source to target
for k in range(n_beams):
    xs, ys, zs = [xs_a[k], ys_a[k], zs_a[k]]
    vect_new = aim_vect[k,:]

    # FIND BASIS VECTOR FOR THE GRIDS, GIVEN THE AIMING VECTOR
    e1,e2 = find_basis(-vect_new)
    # CALCULATE POSITION OF BEAMLETS USING THE BASIS VECTOR
    for i in range(ap_grid[1]): #COLUMNS
        for j in range(ap_grid[0]): #ROWS
            centers[k,i,j,0] = xs+e1[0]*(h_ap_dist*i-tot_size[1]*0.5)+e2[0]*(v_ap_dist*j-tot_size[0]*0.5)
            centers[k,i,j,1] = ys+e1[1]*(h_ap_dist*i-tot_size[1]*0.5)+e2[1]*(v_ap_dist*j-tot_size[0]*0.5)
            centers[k,i,j,2] = zs+e1[2]*(h_ap_dist*i-tot_size[1]*0.5)+e2[2]*(v_ap_dist*j-tot_size[0]*0.5)

    centers_new[k,:,:,:] = centers[k,:,:,:]

    #FIND VECTOR NORMAL TO GRID AFTER ROTATION
    p1 = centers_new[k, 0,0 ,:]
    p2 = centers_new[k,-1,0 ,:]
    p3 = centers_new[k,-1,-1,:]
    v_ort1 = p2-p1
    v_ort2 = p3-p1
    v_ort[k,:] = -np.cross(v_ort1, v_ort2)
    m = math.sqrt(sum(i**2 for i in v_ort[k,:]))
    v_ort[k,:] = np.multiply(v_ort[k,:],1/m)
    
# FIND FOCUS POINT OF TWO SOURCES OF THE SAME INJECTOR
# CHECKING THAT THE AIMING VECTORS CROSS AND THE POINT IS GOOD
for k in range(n_beams/2):
    k2=k+12
    sour1 = [xs_a[k], ys_a[k], zs_a[k]]
    m1=np.linalg.norm(sour1)
    sour2 = [xs_a[k2], ys_a[k2], zs_a[k2]]
    m2=np.linalg.norm(sour2)

    v1 = sour1 + m1*v_ort[k ,:]
    v2 = sour2 + m2*v_ort[k2,:]
    
    A = [sour1[0], sour1[2]]
    B = [v1   [0], v1   [2]]
    C = [sour2[0], sour2[2]]
    D = [v2   [0], v2   [2]]
    x_int, z_int = line_intersection((A,B),(C,D))
    y_int = (x_int-sour1[0])*v_ort[k,1]/v_ort[k,0]+sour1[1]
    intersections[k ,:] = x_int,y_int,z_int
    intersections[k2,:] = x_int,y_int,z_int


#FIND VECTOR FROM BEAMLET GG TO AIMING POINT
for k in range(n_beams):
    aim = intersections[k,:]
    for i in range(ap_grid[1]): #COLUMNS
        for j in range(ap_grid[0]): #ROWS
            origin = centers_new[k,i,j,:]
            new_aimv = aim-origin
            new_aimv = new_aimv/np.linalg.norm(new_aimv)
            v_ort_bmlts[k,i,j,:] = new_aimv

if plot_flag==1:
    xlim=[-15000,15000]
    ylim=[-15000,15000]
    zlim=[-7000,7000]
    #PLOTS
    ##PLOT XZ
    fig2 = plt.figure()
    for i in range(n_beams):
        xs, ys, zs = [xs_a[i], ys_a[i], zs_a[i]]
        plt.plot(centers_new[i,:,:,0], centers_new[i,:,:,2], 'o')
        m=np.linalg.norm([xs,ys,zs])
        plt.plot([xs, xs+v_ort[i,0]*m], [zs, zs+v_ort[i,2]*m],'-')
        if i in (0,12):
            for ii in range(ap_grid[1]): #COLUMNS
                for jj in range(ap_grid[0]): #ROWS
                    x_tmp=centers_new[i,ii,jj,0]
                    z_tmp = centers_new[i,ii,jj,2]
                    plt.plot([x_tmp, x_tmp+v_ort_bmlts[i,ii,jj,0]*m], [z_tmp, z_tmp+v_ort_bmlts[i,ii,jj,2]*m],'-')
                
        plt.text(centers_new[i,12,21,0]+100,centers_new[i,12,21,2]+100,names_label[i],fontsize=15)
        if i<n_beams/2:
          plt.plot(intersections[i,0],intersections[i,2],'x', color='r')

    plt.plot(-R_vess, Z_vess, 'm')
    plt.plot(R_vess, Z_vess, 'm')
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Z')


    ##PLOT XY
    fig4,ax = plt.subplots(1)
    for i in range(n_beams):
        xs, ys, zs = [xs_a[i], ys_a[i], zs_a[i]]
        plt.plot(centers_new[i,:,:,0], centers_new[i,:,:,1], 'o')
        m=np.linalg.norm([xs,ys,zs])
        plt.plot([xs, xs+v_ort[i,0]*m], [ys, ys+v_ort[i,1]*m],'-')
        if i in (0,12):
            for ii in range(ap_grid[1]): #COLUMNS
                for jj in range(ap_grid[0]): #ROWS
                    x_tmp=centers_new[i,ii,jj,0]
                    y_tmp = centers_new[i,ii,jj,1]
                    plt.plot([x_tmp, x_tmp+v_ort_bmlts[i,ii,jj,0]*m], [y_tmp, y_tmp+v_ort_bmlts[i,ii,jj,1]*m],'-')
        if i<n_beams/2:
          plt.plot(intersections[i,0],intersections[i,1],'x',color='r')

    plt.plot([x1_TF,x2_TF,x3_TF,x4_TF,x1_TF],[y1_TF,y2_TF,y3_TF,y4_TF,y1_TF], 'r-')
    plt.text(x3_TF,y3_TF, 'TF COIL 1',fontsize=15, color='red')

    the=np.linspace(0,6.28,100)
    plt.plot(np.min(R_vess)*np.cos(the), np.min(R_vess)*np.sin(the),'m')
    plt.plot(np.max(R_vess)*np.cos(the), np.max(R_vess)*np.sin(the),'m')
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')


    ## ##3D plot
    ## #TOKAMAK
    phi = np.arange(0,2*np.pi,0.02*np.pi)
    x_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
    y_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
    for i,R in enumerate(R_vess):
        x_tok[:,i] = R*np.cos(phi)
        y_tok[:,i] = R*np.sin(phi)
    z_tok = Z_vess

    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(centers_new[:,0,0,0], centers_new[:,0,0,1], centers_new[:,0,0,2], marker='.')
    ax.scatter(centers_new[:,0,-1,0], centers_new[:,0,-1,1], centers_new[:,0,-1,2], marker='.')
    ax.scatter(centers_new[:,-1,0,0], centers_new[:,-1,0,1], centers_new[:,-1,0,2], marker='.')
    ax.scatter(centers_new[:,-1,-1,0], centers_new[:,-1,-1,1], centers_new[:,-1,-1,2], marker='.')

    for k in range(n_beams):
        xs, ys, zs = [xs_a[k], ys_a[k], zs_a[k]]
        mod=np.linalg.norm([xs,ys,zs])
        v_ort_t = v_ort[k,:]*mod
        ax.plot([xs, xs+v_ort_t[0]], [ys, ys+v_ort_t[1]], zs=[zs, zs+v_ort_t[2]], color='red')
        if k<n_beams/2:
            ax.scatter(intersections[k,0], intersections[k,1], intersections[k,2],marker='x')

    ax.plot_surface(x_tok,y_tok,z_tok,color='k',alpha=0.2)

# CALCULATION OF PHI & THETA FOR BBNBI STANDALONE
# WHERE X,Y,Z are the cartesian coordinates in mm
# phi is the angle for the tangency radius (if 0 the injection is normal) (see pietro's figure)
# theta is the vertical tilt
cent_phitheta = np.zeros((n_beams, ap_grid[1], ap_grid[0], 2))
tan_rad_PNB = np.zeros((n_beams))
st_angle = np.zeros((n_beams))
for ii in range(n_beams):
    for jj in range(ap_grid[1]):
       	for kk in range(ap_grid[0]):
            # x1,y1,z1 are the COORDS OF THE BEAMLET ORIGIN
            x1, y1, z1 = centers_new[ii,jj,kk,0], centers_new[ii,jj,kk,1], \
                         centers_new[ii,jj,kk,2]
            #here you must put the angle wrt center
            m1 = y1/x1
            m_aim    = v_ort_bmlts[ii,jj,kk,1]/v_ort_bmlts[ii,jj,kk,0]
            # math.atan2(y, x)
            # Return atan(y / x), in radians. The result is between -pi and pi. 
            # The vector in the plane from the origin to point (x, y) makes this angle with the positive X axis.
            # The point of atan2() is that the signs of both inputs are known to it, so it can compute the correct quadrant for the angle.
            # For example, atan(1) and atan2(1, 1) are both pi/4, but atan2(-1, -1) is -3*pi/4.
            ang1 = math.atan2(y1,x1)
            if ang1<0:
                ang1 = -ang1
            angaim = math.atan2(v_ort_bmlts[ii,jj,kk,1], v_ort_bmlts[ii,jj,kk,0])
            tan_phi = abs((m1-m_aim)/(1+m1*m_aim))
            val=math.atan(tan_phi)

            if names_label[ii] in ('3A','3B','4A','4B','7A','7B','8A','8B'):
                val=-val
            cent_phitheta[ii,jj,kk,0] = val

            z_vec=v_ort_bmlts[ii,jj,kk,2]
            theta = -math.atan2(z_vec, math.sqrt(1-z_vec*z_vec))
            cent_phitheta[ii,jj,kk,1] = theta                

    
# DATA STORAGE FOR BBNBI STANDALONE
# X,Y,Z,PHI,THETA [MM,MM,MM,RAD,RAD]
out_fname = 'PNB_BBNBI.dat'
out_csvname = 'PNB_BBNBI.csv'
outf = open(out_fname,'w')
outcsv = open(out_csvname,'w')
#outf.write("%10s %10s %10s %10s %10s \n" %('X [mm]','Y [mm]','Z [mm]','PHI [rad]','THETA [rad]'))
writer = csv.writer(outcsv, delimiter=',')

for ii in range(n_beams):
    for jj in range(ap_grid[1]):
        for kk in range(ap_grid[0]):
            if jj in (0, ap_grid[1]-1):
                if kk in (0,1,ap_grid[0]-2, ap_grid[0]-1):
                    continue
            elif jj in (1, ap_grid[1]-2):
                if kk in (0, ap_grid[0]-1):
                    continue

                
            x,y,z = centers_new[ii,jj,kk,0], centers_new[ii,jj,kk,1], \
                    centers_new[ii,jj,kk,2]
            phi, theta = cent_phitheta[ii,jj,kk,0], cent_phitheta[ii,jj,kk,1]
            row = [x, y, z, phi, theta]
            writer.writerow(row)
            outf.write(" %10.4f %10.4f %10.4f %10.8f %10.8f \n" %(x, y, z, phi, theta))

outf.close()
outcsv.close()
print "Stored for BBNBI standalone in ", out_fname, " & ", out_csvname

# DATA STORAGE FOR ITM
# X,Y,Z,PHI,THETA [MM,MM,MM,RAD,RAD]
#out_fname = 'PNB_beamlets.dat'
out_csvnameA = 'PNB_A_ITM.csv'
out_csvnameB = 'PNB_B_ITM.csv'
outcsvA = open(out_csvnameA,'w')
outcsvB = open(out_csvnameB,'w')
#outf.write("%10s %10s %10s %10s %10s \n" %('X [mm]','Y [mm]','Z [mm]','PHI [rad]','THETA [rad]'))
#writerA = csv.writer(outcsvA, delimiter=' ')
#writerB = csv.writer(outcsvB, delimiter=' ')

for ii in range(n_beams):
    for jj in range(ap_grid[1]):
        for kk in range(ap_grid[0]):
            if jj in (0, ap_grid[1]-1):
                if kk in (0,1,ap_grid[0]-2, ap_grid[0]-1):
                    continue
            elif jj in (1, ap_grid[1]-2):
                if kk in (0, ap_grid[0]-1):
                    continue

                
            x,y,z = centers_new[ii,jj,kk,0], centers_new[ii,jj,kk,1], \
                    centers_new[ii,jj,kk,2]
            vx,vy,vz = v_ort_bmlts[ii,jj,kk,0], v_ort_bmlts[ii,jj,kk,1],\
                       v_ort_bmlts[ii,jj,kk,2]

            phi = math.atan2(vy,vx)
            theta = cent_phitheta[ii,jj,kk,1]
            
            row = np.array((x*1e-3, y*1e-3, z*1e-3, theta, phi))
            
            if ii < n_beams/2:
                np.savetxt(outcsvA, row.reshape(1,row.shape[0]), fmt='%#010.6f', delimiter=' ')
            else:
                np.savetxt(outcsvB, row.reshape(1,row.shape[0]), fmt='%#010.6f', delimiter=' ')

outcsvA.close()
outcsvB.close()

print "Stored for ITM in ", out_csvnameA , " & ", out_csvnameB





plt.show()
