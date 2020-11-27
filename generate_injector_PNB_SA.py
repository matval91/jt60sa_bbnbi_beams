import math
import os

inj=1
npinis=4
nbeamlets=1020

energy=85e3
frac1=0.78
frac2=0.15
frac3=0.07
isotope=2
power=1e6
div=0.01754

output_R=4.5
output_phi=130.0/180.0*math.pi
output_w=2.0
output_h=8.0

if not os.path.exists("inj%d" % inj):
  os.makedirs("inj%d" % inj)

beamlets=[]
f=open("beamlets.csv","r")
lines=f.readlines()
for line in lines:
  parts=line.split(",")
  beamlets.append([float(i) for i in parts])
f.close()

f=open("inj%d/injector" % inj,"w")
f.write("""\
# Injector [vMarch2010] (Pinis, obstacles, output surface, Pini powers):
       %d          # Injector: # Injector ID
       %d          # Injector: # of Pinis
""" % (inj, npinis))
k=0
for i in range(1,npinis+1):
  f.write("""\
# Pini:
  %d2           # Beamlet weights id
  %d1           # Energy fractions id
 0.00000000000000D+00   0.00000000000000D+00      # Horizontal & vertical misalignment
      %d             # Pini: # of beamlets
""" % (inj, inj, nbeamlets))
  for j in range(0,nbeamlets):
    f.write("""\
# Beamlet (IDs, coords, direction):
  %d5    %d3 # Beamlet: disp. id, ANums id.
 %.14E %.14E %.14E  # Point: x, y, z
 %.14E %.14E   0.10000000000000D+01  # Vector: theta, phi, length
""" % (inj, inj, beamlets[k][0]/1e3, beamlets[k][1]/1e3, beamlets[k][2]/1e3, beamlets[k][4]+math.pi/2, beamlets[k][3]+math.pi+math.atan2(beamlets[k][1],beamlets[k][0])))
    k=k+1

phi_min=output_phi-math.atan(0.5*output_w/output_R)
phi_max=output_phi+math.atan(0.5*output_w/output_R)
x1=output_R*math.cos(phi_min)
x2=output_R*math.cos(phi_max)
y1=output_R*math.sin(phi_min)
y2=output_R*math.sin(phi_max)
z1=-0.5*output_h
z2=0.5*output_h

f.write("""\
       0           # Injector: # of obstacles
# Output surface:
# Quad (2 Triangles):
# Triangle (3 Points):
%.14E %.14E %.14E # Point: x, y, z
%.14E %.14E %.14E # Point: x, y, z
%.14E %.14E %.14E # Point: x, y, z
# Triangle (3 Points):
%.14E %.14E %.14E # Point: x, y, z
%.14E %.14E %.14E # Point: x, y, z
%.14E %.14E %.14E # Point: x, y, z
  %d4           # Pini powers id
""" % (x1, y1, z1, x2, y2, z1, x1, y1, z2, x1, y1, z2, x2, y2, z1, x2, y2, z2, inj) )
f.close()

f=open("inj%d/fullenergies" % inj, "w")
f.write("# Injector: Full energies of particles from each PINI\n")
for i in range(1,npinis+1):
  f.write("%.14E # PINI %d\n" % (energy, i))
f.close()

f=open("inj%d/distributions" % inj, "w")
f.write("          4  - # of discrete distributions\n")
f.write("""\
# The possible energy fractions & probabilities                                 
  %d1       3              # DiscreteDistribution: id, # of keys
 0.10000000000000D+01   %.14E      # DiscreteDistribution: data point
 0.50000000000000D+00   %.14E      # DiscreteDistribution: data point
 0.33333333333333D+00   %.14E      # DiscreteDistribution: data point
""" % (inj, frac1, frac2, frac3))
f.write("""\
# Relative weights of the beamlets in each PINI. 
  %d2     %d              # DiscreteDistribution: id, # of keys
""" % (inj, nbeamlets))
for i in range(1,nbeamlets+1):
  f.write("%d 1      # DiscreteDistribution: data point\n" % i)
f.write("""\
# The possible isotopes & probabilities of emitted particles                    
  %d3       1              # DiscreteDistribution: id, # of keys
 %.14E   0.10000000000000D+01      # DiscreteDistribution: data point
""" % (inj, isotope))
f.write("""\
# The powers of the PINIs in Watts
  %d4       %d              # DiscreteDistribution: id, # of keys
""" % (inj, npinis))
for i in range(1,npinis+1):
  f.write("%.14E %.14E      # DiscreteDistribution: data point\n" % (i, power))
f.write("""\
           1  - # of distributions
# Beamlet dispersion distribution:
%d5  # Distribution id
140  # Distribution: # of data points
""" % inj)
for i in range(0, 140):
  x=i/140.0*0.06
  f.write("%.14E %.14E  # Distribution data point\n" % (x, ((1.0/(math.pi*div**2))*math.exp(-((x/div)**2)))*2*math.pi*math.sin(x)))
f.close()
