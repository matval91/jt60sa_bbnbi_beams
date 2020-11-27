# Script to generate the beamlet data for the JT60-SA NEGATIVE NBs
# with data obtained from the references in the JT60-SA document
# M. Vallar - 11/2016 matteo.vallar@igi.cnr.it
#-----------------------------------------------------------------

# MODULES TO LOAD
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
###################################################
plot_flag=1

def rotate_point(x,y,theta):
    """
    Calculate the rotation of a point given the angle theta
    """
    xnew=x*math.cos(theta)-y*math.sin(theta)
    ynew=x*math.sin(theta)+y*math.cos(theta)
    return xnew, ynew 

def find_basis(vector):
    """
    Find the basis for the plane defined by vector. One vector is in the form (a,b,0) (2D vector in 3D space), with a*a+b*b=1
    """
    b=vector[0]/math.sqrt(vector[0]*vector[0]+vector[1]*vector[1])
    a=-b*vector[1]/vector[0]
    basis_1=[a,b,0]
    basis_2=np.cross(basis_1,vector)
    basis_2=basis_2/np.linalg.norm(basis_2)
    return basis_1, basis_2


# DESIGN PARAMETERS
n_seg = 5 #Number of segments (i.e. grounded grids) in each N-NB injector.

# DATA RELATIVE TO THE SEGMENTS
seg_size = [168, 437]   # segment size in mm [row,col]
ap_grid = [9,24] # number of apertures in segment [row, col]
h_ap_dist = 19 #horizontal distance bwt two beamlets center
v_ap_dist = 21 # vertical  ""    ""   ""   ""   ""   ""
#-------------------------------------------

# DATA RELATIVE TO THE FULL INJECTOR
R_foc = float(23660)
#-----------------------------------
#VESSEL coordinates (just for plotting)
vess_fname='2D_vess.txt'
vess=np.loadtxt(vess_fname, dtype=float, unpack=True)
R_vess, Z_vess = vess[0,:]*1e3, vess[1,:]*1e3
#COORDINATES OF TF COIL1: (for plotting, in order to check where is the correct 0 of toroidal angle)
# this is as if the TF coil1 is at +x axis
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
#-----------------------------------

#COORDINATES OF CENTERS OF SOURCES Up AND Low
xs =  27751.20
ys = -2855.287
Z_U = 590.5
Z_L = -1690.5

#ROTATION THAT ACCOUNTS FOR THE CORRECT TOROIDAL 0
r_source = math.sqrt(xs*xs+ys*ys)
original_angle = math.atan2(ys,xs)
xs = r_source*math.cos(original_angle+ang_rotation)
ys = r_source*math.sin(original_angle+ang_rotation)
#=======================

#From excel file
ang_xy_U = math.pi-0.5*math.pi/180.+ang_rotation #angle on xy plane
ang_xy_L = ang_xy_U
ang_xz_U = -2.75*math.pi/180. #angle on rz plane
ang_xz_L = 2.75*math.pi/180.

#FIND AIMING VECTORS FOR EACH SEGMENT:
alpha_g=np.array([1,0.5,0,-0.5,-1]) # From word file
alpha_r = np.multiply(alpha_g,math.pi/180.) #conversion to rad
vec_aim_U=np.zeros((5,3),dtype=float) #array of aiming vectors UP
vec_aim_L=np.zeros((5,3),dtype=float) #array of aiming vectors LOW
for ii in range(5):
    x = math.cos(ang_xy_U)*math.cos(ang_xz_U+alpha_r[ii])
    y = math.sin(ang_xy_U)*math.cos(ang_xz_U+alpha_r[ii])
    z = math.sin(ang_xz_U+alpha_r[ii])
    vec_aim_U[ii,:]=np.array((x,y,z))
    x = math.cos(ang_xy_L)*math.cos(ang_xz_L+alpha_r[ii])
    y = math.sin(ang_xy_L)*math.cos(ang_xz_L+alpha_r[ii])
    z = math.sin(ang_xz_L+alpha_r[ii])
    vec_aim_L[ii,:]=np.array((x,y,z))


#FINDING POSITION OF GRID CENTERS FOR EACH SEGMENT
gridcent_U = np.zeros((5,3), dtype=float) # array of the centers of the segments UP
gridcent_L = np.zeros((5,3), dtype=float) # array of the centers of the segments LOW
R=R_foc
for ii in range(5):
    angle_source=math.atan2(ys,xs)
    dR=R*(1-math.cos(alpha_r[ii]))
    dx=dR*math.cos(ang_xy_U)
    dy=dR*math.sin(ang_xy_U)
    gridcent_U[ii,0]=xs+dx
    gridcent_U[ii,1]=ys+dy
    gridcent_U[ii,2]=Z_U+float(230)*(ii-2)
    gridcent_L[ii,0]=xs+dx
    gridcent_L[ii,1]=ys+dy
    gridcent_L[ii,2]=Z_L+float(230)*(ii-2)



#CALCULATIONS OF CENTERS OF THE BEAMLETS
centers = np.zeros((ap_grid[1], ap_grid[0], 3), dtype=float) #temporary variable
centers_U = np.zeros((5,ap_grid[1], ap_grid[0], 3), dtype=float) #beamlets x,y,z position
centers_L = np.zeros((5,ap_grid[1], ap_grid[0], 3), dtype=float) #beamlets x,y,z position

# VECTORS NORMAL TO THE GRID
v_ort_U = np.zeros((5,3), dtype=float)
v_ort_L = np.zeros((5,3), dtype=float)
for k in range(n_seg):
   for pos in ('u','l'):
        if pos=='u':
            vect_new = vec_aim_U[k]
            xs, ys, zs = [gridcent_U[k,0], gridcent_U[k,1], gridcent_U[k,2]]
        else:
            vect_new = vec_aim_L[k]
            xs, ys, zs = [gridcent_L[k,0], gridcent_L[k,1], gridcent_L[k,2]]      
        # FIND BASIS VECTOR FOR THE GRIDS, GIVEN THE AIMING VECTOR
        e1,e2 = find_basis(-vect_new)
        # CALCULATE POSITION OF BEAMLETS USING THE BASIS VECTOR
        for i in range(ap_grid[1]): #COLUMNS
            for j in range(ap_grid[0]): #ROWS
                centers[i,j,0] = xs+e1[0]*(h_ap_dist*i-seg_size[1]*0.5)+e2[0]*(v_ap_dist*j-seg_size[0]*0.5)
                centers[i,j,1] = ys+e1[1]*(h_ap_dist*i-seg_size[1]*0.5)+e2[1]*(v_ap_dist*j-seg_size[0]*0.5)
                centers[i,j,2] = zs+e1[2]*(h_ap_dist*i-seg_size[1]*0.5)+e2[2]*(v_ap_dist*j-seg_size[0]*0.5)

        if pos=='u':
            centers_U[k,:,:,:] = centers[:,:,:]
        else:
            centers_L[k,:,:,:] = centers[:,:,:]
        centers = np.zeros((ap_grid[1], ap_grid[0], 3), dtype=float) #re-initialising to 0

   #FIND VECTOR NORMAL TO GRID AFTER ROTATION
   p1 = centers_U[k, 0,0 ,:]
   p2 = centers_U[k,-1,0 ,:]
   p3 = centers_U[k,-1,-1,:]
   v_ort1 = p2-p1
   v_ort2 = p3-p1
   v_ort_U[k,:] = np.cross(v_ort1, v_ort2)
   m=np.linalg.norm(v_ort_U[k,:])
   v_ort_U[k,:] = np.multiply(v_ort_U[k,:],1/m)
   
   p1 = centers_L[k, 0,0 ,:]
   p2 = centers_L[k,-1,0 ,:]
   p3 = centers_L[k,-1,-1,:]
   v_ort1 = p2-p1
   v_ort2 = p3-p1
   v_ort_L[k,:] = np.cross(v_ort1, v_ort2)
   m=np.linalg.norm(v_ort_L[k,:])
   v_ort_L[k,:] = np.multiply(v_ort_L[k,:],1/m)

##PLOTS
if plot_flag==1:
    #PLOT XZ
    fig2 = plt.figure()
    for i in range(n_seg):
        for pos in ('u','l'):
            mod = R_foc
            if pos == 'u':
                xst,yst,zst = gridcent_U[i,:]
                plt.plot(centers_U[i,:,:,0], centers_U[i,:,:,2], '.')
                plt.plot([xst, xst+v_ort_U[i,0]*mod], [zst, zst+v_ort_U[i,2]*mod],'-')
            if pos == 'l':
                xst,yst,zst = gridcent_L[i,:]
                plt.plot(centers_L[i,:,:,0], centers_L[i,:,:,2], '.')
                plt.plot([xst, xst+v_ort_L[i,0]*mod], [zst, zst+v_ort_L[i,2]*mod],'-')

    plt.plot(R_vess, Z_vess,'m')
    plt.plot(-R_vess,Z_vess,'m')
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Z')

    #PLOT XY
    fig2, ax = plt.subplots(1)
    for i in range(n_seg):
        for pos in ('u','l'):
            mod = R_foc
            if pos == 'u':
                xst,yst,zst = gridcent_U[i,:]
                plt.plot(centers_U[i,:,:,0], centers_U[i,:,:,1], '.')
                plt.plot([xst, xst+v_ort_U[i,0]*mod], [yst, yst+v_ort_U[i,1]*mod],'-')
            if pos == 'l':
                xst,yst,zst = gridcent_L[i,:]
                plt.plot(centers_L[i,:,:,0], centers_L[i,:,:,1], '.')
                plt.plot([xst, xst+v_ort_L[i,0]*mod], [yst, yst+v_ort_L[i,1]*mod],'-')
    plt.plot([x1_TF,x2_TF,x3_TF,x4_TF,x1_TF],[y1_TF,y2_TF,y3_TF,y4_TF,y1_TF], 'r-')
    plt.text(x3_TF,y3_TF, 'TF COIL 1',fontsize=15, color='red')

    the=np.linspace(0,6.28,100)
    plt.plot(np.min(R_vess)*np.cos(the), np.min(R_vess)*np.sin(the),'m')
    plt.plot(np.max(R_vess)*np.cos(the), np.max(R_vess)*np.sin(the),'m')
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')

    ##3D plot
    #TOKAMAK
    #shape of tokamak
    phi = np.arange(0,2*np.pi,0.02*np.pi)
    x_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
    y_tok = np.zeros((len(phi),len(R_vess)),dtype=float)
    for i,R in enumerate(R_vess):
        x_tok[:,i] = R*np.cos(phi)
        y_tok[:,i] = R*np.sin(phi)
    z_tok = Z_vess


    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for k in range(n_seg):
        for pos in ('u','l'):
            m=R_foc
            if pos == 'u':
                xs, ys, zs = gridcent_U[k,:]
                ax.plot([xs, xs+v_ort_U[k, 0]*m], [ys, ys+v_ort_U[k, 1]*m], zs=[zs, zs+v_ort_U[k, 2]*m], color='red')
            if pos == 'l':
                xs, ys, zs = gridcent_L[k,:]
                ax.plot([xs, xs+v_ort_L[k, 0]*m], [ys, ys+v_ort_L[k, 1]*m], zs=[zs, zs+v_ort_L[k, 2]*m], color='red')

    ax.plot_surface(x_tok,y_tok,z_tok,color='k',alpha=0.2)
    ax.axis('equal')


# EXPORT DATA TO FILE
# CALCULATIONS FOR HAVING THE OUTPUT IN X,Y,Z,PHI,THETA
# WHERE X,Y,Z are the cartesian coordinates in mM
# phi is the angle for the tangency radius (if 0 the injection is normal) (see pietro's figure)
# theta is the vertical tilt
cent_phitheta_U = np.zeros((n_seg, ap_grid[1], ap_grid[0], 2))
cent_phitheta_L = np.zeros((n_seg, ap_grid[1], ap_grid[0], 2))

# AIMING ANGLE ON XY OF U SOURCE FROM DATA
x1 = centers_U[4, 12, 4, 0]
y1 = centers_U[4, 12, 4, 1]
x2 = centers_U[4, -1, 4, 0]
y2 = centers_U[4, -1, 4, 1]
xprime = centers_U[4,11,4,0]
yprime = centers_U[4,11,4,1]
m=float(y2-y1)/(x2-x1)
m_aim_U = -1/m

# AIMING ANGLE ON XY  OF L SOURCE FROM DATA
x1 = centers_L[4, 12, 4, 0]
y1 = centers_L[4, 12, 4, 1]
x2 = centers_L[4, -1, 4, 0]
y2 = centers_L[4, -1, 4, 1]
xprime = centers_U[4,11,4,0]
yprime = centers_U[4,11,4,1]
m=float(y2-y1)/(x2-x1)
m_aim_L = -1/m

tan_rad_NNB = []
# CALCULATION OF PHI & THETA
for ii in range(n_seg):
	for jj in range(ap_grid[1]):
		for kk in range(ap_grid[0]):
			x1, y1, z1 = centers_U[ii,jj,kk,0], centers_U[ii,jj,kk,1], \
                                     centers_U[ii,jj,kk,2]
			m1 = y1/x1
			tan_phi = abs((m1-m_aim_U)/(1+m1*m_aim_U))
			cent_phitheta_U[ii,jj,kk,0] = math.atan(tan_phi)
			theta = -math.atan2(v_ort_U[ii,2],math.sqrt(1-v_ort_U[ii,2]*v_ort_U[ii,2]))
			cent_phitheta_U[ii,jj,kk,1] = theta			

			x1, y1, z1 = centers_L[ii,jj,kk,0], centers_L[ii,jj,kk,1], \
                                     centers_L[ii,jj,kk,2]
			m1 = y1/x1
			tan_phi = abs((m1-m_aim_L)/(1+m1*m_aim_L))
			cent_phitheta_L[ii,jj,kk,0] = math.atan(tan_phi)

			theta = -math.atan2(v_ort_L[ii,2],math.sqrt(1-v_ort_L[ii,2]*v_ort_L[ii,2]))
			cent_phitheta_L[ii,jj,kk,1] = theta	

# DATA STORAGE
# X,Y,Z,PHI,THETA [MM,MM,MM,RAD,RAD]
out_fname = 'NNB_BBNBI.dat'
out_csvname = 'NNB_BBNBI.csv'
outf = open(out_fname,'w')
outcsv = open(out_csvname,'w')
writer = csv.writer(outcsv, delimiter=',')
for pos in ('u','l'):
	for ii in range(n_seg):
		for jj in range(ap_grid[1]):
			for kk in range(ap_grid[0]):
				if pos=='u':
					x,y,z = centers_U[ii,jj,kk,0], centers_U[ii,jj,kk,1], \
			    	 		centers_U[ii,jj,kk,2]
					phi, theta = cent_phitheta_U[ii,jj,kk,0], cent_phitheta_U[ii,jj,kk,1]
				else:
					x,y,z = centers_L[ii,jj,kk,0], centers_L[ii,jj,kk,1], \
			    	 		centers_L[ii,jj,kk,2]
					phi, theta = cent_phitheta_L[ii,jj,kk,0], cent_phitheta_L[ii,jj,kk,1]
				outf.write(" %10.4f %10.4f %10.4f %10.4f %10.4f \n" %(x, y, z, phi, theta))
                                row = [x, y, z, phi, theta]
                                writer.writerow(row)
outf.close()
outcsv.close()
print "NNBS beamlets for BBNBI standalone are stored in ", out_fname, " & ", out_csvname



# DATA STORAGE FOR ITM
# X,Y,Z,PHI,THETA [MM,MM,MM,RAD,RAD]
out_csvname = 'NNB_ITM.csv'
outcsv = open(out_csvname,'wab')
for pos in ('u','l'):
	for ii in range(n_seg):
		for jj in range(ap_grid[1]):
			for kk in range(ap_grid[0]):
				if pos=='u':
					x,y,z = centers_U[ii,jj,kk,0], centers_U[ii,jj,kk,1], \
			    	 		centers_U[ii,jj,kk,2]
                                        phi   = math.atan2(v_ort_U[ii,1],v_ort_U[ii,0])
                                        theta = math.atan2(math.sqrt(1-v_ort_U[ii,2]*v_ort_U[ii,2]),v_ort_U[ii,2])
				else:
					x,y,z = centers_L[ii,jj,kk,0], centers_L[ii,jj,kk,1], \
			    	 		centers_L[ii,jj,kk,2]
                                        phi   = math.atan2(v_ort_L[ii,1],v_ort_L[ii,0])
                                        theta = math.atan2(math.sqrt(1-v_ort_L[ii,2]*v_ort_L[ii,2]),v_ort_L[ii,2])
                                row = np.array((x*1e-3, y*1e-3, z*1e-3, theta, phi))
                                #writer.writerow(" %10.4f %10.4f %10.4f %10.4f %10.4f \n" %(x, y, z, phi, theta))
                                np.savetxt(outcsv, row.reshape(1,row.shape[0]), fmt='%#010.6f', delimiter=' ')
                                

outcsv.close()
print "NNBS beamlets for ITM are stored in ", out_csvname
plt.show()
