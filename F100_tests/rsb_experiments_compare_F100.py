import numpy as np
from netCDF4 import Dataset
import seawater as sw
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import depths as depth
import pandas as pd
from matplotlib.ticker import LogLocator,Locator
import glob as glob

# 1m vs 3m vs 10 m

res0 = 1
res1 = 3 # horizontal resolution in m
res2 = 10
res4 = 10
res5 = 10
res6 = 3

cent0 = 143 # pipe center + grid points to edge of pipe
#cent1 = 261
cent2 = 64
cent1 = 64
cent4 = 64
cent5 = 64
cent6 = 261

l_cen0 = 256
l_cen1 = 256
l_cen2 = 64
l_cen4 = 64
l_cen5 = 64
l_cen6 = 256



tpas = 1500

w_p0 = 30/2 # width of half of line source in grid points
w_p1 = 10/2
w_p2 = 2 # 3 grid points wide
w_p3 = 0 # 5 grid points wide
w_p4 = 2
w_p5 = 2
w_p6 = 10/2

l_p0 = 500/2
l_p1 = 300/2  # length of half of pipe in grid points
l_p2 = 90/2
l_p4 = 90/2
l_p5 = 90/2
l_p6 = 300/2

# water column depth 
w_col = 60.

# effluent properties
temp_e = 26.9
salt_e = 1.2

'''
# netCDF outputs --> change how this works to grab diff files from folders (glob)
files0 = list(sorted(glob.glob('/maya/minnaho/croco_1m_F100/*')))
files1 = list(sorted(glob.glob('/maya/minnaho/croco_10m_F100/*')))
files2 = list(sorted(glob.glob('/maya/minnaho/croco_F100_dt.5/*')))
files3 = list(sorted(glob.glob('/maya/minnaho/croco_F100_5x300/*')))
#files4 = list(sorted(glob.glob('/maya/minnaho/croco_F100_dt.5/*')))    

# first files with grid variables and first few time steps
output0_0 = files0[0] 
output0_1 = files1[0] 
output0_2 = files2[0]
output0_3 = files3[0]
'''
output0_4 = 'pipe_10m_nh_F100_gray_max_grd.nc'
output0_3 = 'pipe_10m_nh_F100_kpp_grd.nc'
output0_5 = 'pipe_10m_nh_F100_noslip_grd.nc'
output0_6 = 'pipe_3m_dt1_nh_F100_gray_grd.nc'
# last file or second to last for time averaging
'''
output_0 = files0[-2]
output_1 = files1[-1]
output_2 = files2[-1]
output_3 = files3[-1]
'''
output_4 = 'pipe_10m_nh_F100_gray_max_his.20000101000000.nc'
output_3 = 'pipe_10m_nh_F100_kpp_his.20000101000000.nc'
output_5 = 'pipe_10m_nh_F100_noslip_his.20000101000000.nc'
output_6 = 'pipe_3m_dt1_nh_F100_gray_his.20000101000000.nc'
#output_4 = files4[-1]

# time to average over
stime = 10
#etime = 24
etime = 21 # ran 3m gray mix F100 until 21 outputs

# for plotting linear-logarithmic x axis
class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in xrange(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))

###########################
# find Froude # F = (u^3)/b
###########################
# temp and salinity profiles
temp_nc = np.array(Dataset(output_4,'r').variables['temp'])
salt_nc = np.array(Dataset(output_4,'r').variables['salt'])
dens_a = sw.dens0(salt_nc[0,0,cent1,cent1],temp_nc[0,0,cent1,cent1]) # density of ambient
                                                                 # at pipe (bottom)
dens_e = sw.dens0(salt_e,temp_e) # density of effluent

g = 9.8
Q = 10. # m3/s
L = 900. # m
q = Q/L
buoy_roms = (g*(dens_a-dens_e)/dens_a)*q
'''
# get u values
u_nc0 = np.array(Dataset(output0_0,'r').variables['u'][0,0,0,0])
u_nc1 = np.array(Dataset(output0_1,'r').variables['u'][0,0,0,0])
u_nc2 = np.array(Dataset(output0_2,'r').variables['u'][0,0,0,0])
u_nc3 = np.array(Dataset(output0_3,'r').variables['u'][0,0,0,0])
#u_nc4 = np.array(Dataset(output0_4,'r').variables['u'][0,0,0,0])

u_arr = np.array([u_nc0,u_nc1,u_nc2,u_nc3])
'''
#F_plt = ((u_arr)**3)/buoy_roms
F_plt = np.array([100,100,100,100])

#######################################
# find height of maximum concentration
# at different distances away from pipe
#######################################
#tpas_nc0 = np.array(Dataset(output_0,'r').variables['tpas'])
#tpas_nc1 = np.array(Dataset(output_1,'r').variables['tpas'])
#tpas_nc2 = np.array(Dataset(output_2,'r').variables['tpas'])
tpas_nc3 = np.array(Dataset(output_3,'r').variables['trace1'])
tpas_nc4 = np.array(Dataset(output_4,'r').variables['trace1'])
tpas_nc5 = np.array(Dataset(output_5,'r').variables['trace1'])
tpas_nc6 = np.array(Dataset(output_6,'r').variables['trace1'])

# linearly stratified so N is constant
z0 = 60
z1 = z0
z2 = z0
z3 = z0
z4 = z0
z5 = z0
z6 = z0
dens_0 = sw.dens0(salt_nc[0,-1,cent1,cent1],temp_nc[0,-1,cent1,cent1]) # density at surface

N0 = ((-g/dens_a)*((dens_0-dens_a)/z0))**.5
N1 = ((-g/dens_a)*((dens_0-dens_a)/z1))**.5
N2 = ((-g/dens_a)*((dens_0-dens_a)/z2))**.5
N3 = ((-g/dens_a)*((dens_0-dens_a)/z3))**.5
N4 = ((-g/dens_a)*((dens_0-dens_a)/z4))**.5
N5 = ((-g/dens_a)*((dens_0-dens_a)/z5))**.5
N6 = ((-g/dens_a)*((dens_0-dens_a)/z6))**.5
#N4 = N3

# lb should be same for all because constant N
lb0 = (buoy_roms**(1/3.))/N0
lb1 = (buoy_roms**(1/3.))/N1
lb2 = (buoy_roms**(1/3.))/N2
lb3 = (buoy_roms**(1/3.))/N3
lb4 = (buoy_roms**(1/3.))/N4
lb5 = (buoy_roms**(1/3.))/N5
lb6 = (buoy_roms**(1/3.))/N6

# find top and max concentration of wastefield by 
# taking average profiles of tpas for highest concentration
# at last time step at different points 
# (0, 100, 200, 300, 400, 500, 600, 700 m away from pipe)
# over length of pipe (cent1-l_p:cent1+l_p)
#l_real = np.array([0.,100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500.,1600.,1700.,1800.,1900.,2000.]) # m away from pipe
#l_real = np.array([0.,100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500.]) # m away from pipe
#l_dist = (l_real/res).astype(int) # grid points away from pipe
l_dist0 = np.array([list(range(1,700))])[0]
l_dist1 = np.array([list(range(1,700))])[0]
l_dist2 = np.array([list(range(1,191))])[0]
l_dist3 = np.array([list(range(1,192))])[0]
l_dist4 = np.array([list(range(1,192))])[0]
l_dist5 = np.array([list(range(1,192))])[0]
l_dist6 = np.array([list(range(1,700))])[0]
l_real0 = l_dist0*res0
l_real1 = l_dist1*res1
l_real2 = l_dist2*res2
l_real3 = l_dist3*res2
l_real4 = l_dist4*res4
l_real5 = l_dist5*res2
l_real6 = l_dist6*res6

# after the following, left with vertical profile averaged over length of pipe at last time step
# for each distance away from pipe
#avg_prof0 = np.empty((len(l_dist0),tpas_nc0.shape[1])) # shape (l_dist.shape[0],64)
#avg_prof1 = np.empty((len(l_dist2),tpas_nc1.shape[1])) # shape (l_dist.shape[0],64)
#avg_prof2 = np.empty((len(l_dist1),tpas_nc2.shape[1])) # shape (l_dist.shape[0],64)
avg_prof3 = np.empty((len(l_dist3),tpas_nc3.shape[1])) # shape (l_dist.shape[0],64)
avg_prof4 = np.empty((len(l_dist4),tpas_nc4.shape[1])) # shape (l_dist.shape[0],64)
avg_prof5 = np.empty((len(l_dist5),tpas_nc5.shape[1])) # shape (l_dist.shape[0],65)
avg_prof6 = np.empty((len(l_dist6),tpas_nc6.shape[1])) 

for l_i in range(len(l_dist4)):
    avg_prof4[l_i] = np.mean(np.mean(tpas_nc4[stime:etime,:,int(l_cen4-l_p4):int(l_cen4+l_p4),int(cent4+w_p4+l_dist4[l_i])],axis=0),axis=1)

for l_i in range(len(l_dist3)):
    avg_prof3[l_i] = np.mean(np.mean(tpas_nc3[stime:etime,:,int(l_cen2-l_p2):int(l_cen2+l_p2),int(cent2+w_p2+l_dist3[l_i])],axis=0),axis=1)
    avg_prof5[l_i] = np.mean(np.mean(tpas_nc5[stime:etime,:,int(l_cen2-l_p2):int(l_cen2+l_p2),int(cent2+w_p2+l_dist5[l_i])],axis=0),axis=1)

for l_i in range(len(l_dist6)):
    avg_prof6[l_i] = np.mean(np.mean(tpas_nc6[stime:etime,:,int(l_cen6-l_p6):int(l_cen6+l_p6),int(cent6+w_p6+l_dist6[l_i])],axis=0),axis=1)

'''

for l_i in range(len(l_dist0)):
    avg_prof0[l_i] = np.mean(np.mean(tpas_nc0[stime:etime,:,l_cen0-l_p0:l_cen0+l_p0,cent0+w_p0+l_dist0[l_i]],axis=0),axis=1)

for l_i in range(len(l_dist2)):
    avg_prof1[l_i] = np.mean(np.mean(tpas_nc1[stime:etime,:,l_cen2-l_p2:l_cen2+l_p2,cent2+w_p2+l_dist2[l_i]],axis=0),axis=1)

for l_i in range(len(l_dist1)):
    avg_prof2[l_i] = np.mean(np.mean(tpas_nc2[stime:etime,:,l_cen1-l_p1:l_cen1+l_p1,cent1+w_p1+l_dist1[l_i]],axis=0),axis=1)
    avg_prof3[l_i] = np.mean(np.mean(tpas_nc3[stime:etime,:,l_cen1-l_p1:l_cen1+l_p1,cent1+w_p3+l_dist1[l_i]],axis=0),axis=1)
#    avg_prof4[l_i] = np.mean(np.mean(tpas_nc4[stime:etime,:,l_cen1-l_p1:l_cen1+l_p1,cent1+w_p1+l_dist1[l_i]],axis=0),axis=1)

maxzm0 = np.empty((len(range(l_cen0-l_p0,l_cen0+l_p0)),len(l_dist0)))
maxzm1 = np.empty((len(range(l_cen2-l_p2,l_cen2+l_p2)),len(l_dist2)))
maxzm2 = np.empty((len(range(l_cen1-l_p1,l_cen1+l_p1)),len(l_dist1)))
maxzm3 = np.empty((len(range(l_cen1-l_p1,l_cen1+l_p1)),len(l_dist1)))

# height to minimum dilution (max concentration)
# get depths of s_rho levels at last time step
[z_sigmas0,Cs] = depth.get_depths(output_0,output0_0,-1,'w','new')
[z_sigmas1,Cs] = depth.get_depths(output_1,output0_1,-1,'w','new')
[z_sigmas2,Cs] = depth.get_depths(output_2,output0_2,-1,'w','new')
[z_sigmas3,Cs] = depth.get_depths(output_3,output0_3,-1,'w','new')
[z_sigmas4,Cs] = depth.get_depths(output_4,output0_4,-1,'w','new')

print('zm')
for ldist_i in range(len(l_dist0)):
    for len_i in range(l_cen0-l_p0,l_cen0+l_p0):
        maxzm0[len_i-(l_cen0-l_p0),ldist_i] = z_sigmas0[np.where(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0)==max(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0))),len_i,l_dist0[ldist_i]][0][0]

for ldist_i in range(len(l_dist2)):
    for len_i in range(l_cen2-l_p2,l_cen2+l_p2):
        maxzm1[len_i-(l_cen2-l_p2),ldist_i] = z_sigmas1[np.where(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0)==max(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0))),len_i,l_dist2[ldist_i]][0][0]

print('zm1')
for ldist_i in range(len(l_dist1)):
    for len_i in range(l_cen1-l_p1,l_cen1+l_p1):
        maxzm2[len_i-(l_cen1-l_p1),ldist_i] = z_sigmas2[np.where(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0)==max(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0))),len_i,l_dist1[ldist_i]][0][0]
        maxzm3[len_i-(l_cen1-l_p1),ldist_i] = z_sigmas3[np.where(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0)==max(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0))),len_i,l_dist1[ldist_i]][0][0]

mean_zm = (np.array([np.mean(maxzm0,axis=0),np.mean(maxzm1,axis=0),np.mean(maxzm2,axis=0),np.mean(maxzm3,axis=0)])+w_col)/lb0

# percent of maximum to calculate thickness and top of plume
perc = 0.05
maxhe0 = np.empty((len(range(l_cen0-l_p0,l_cen0+l_p0)),len(l_dist0)))
maxhe1 = np.empty((len(range(l_cen2-l_p2,l_cen2+l_p2)),len(l_dist2)))
maxhe2 = np.empty((len(range(l_cen1-l_p1,l_cen1+l_p1)),len(l_dist1)))
maxhe3 = np.empty((len(range(l_cen1-l_p1,l_cen1+l_p1)),len(l_dist1)))

maxze0 = np.empty((len(range(l_cen0-l_p0,l_cen0+l_p0)),len(l_dist0)))
maxze1 = np.empty((len(range(l_cen2-l_p2,l_cen2+l_p2)),len(l_dist2)))
maxze2 = np.empty((len(range(l_cen1-l_p1,l_cen1+l_p1)),len(l_dist1)))
maxze3 = np.empty((len(range(l_cen1-l_p1,l_cen1+l_p1)),len(l_dist1)))

print('he')
for ldist_i in range(len(l_dist0)):
    for len_i in range(l_cen0-l_p0,l_cen0+l_p0):
        maxhe0[len_i-(l_cen0-l_p0),ldist_i] = z_sigmas0[np.where(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0)>max(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist0[ldist_i]] - z_sigmas0[np.where(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0)>max(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0))*perc)[0][0],len_i,l_dist0[ldist_i]]
        maxze0[len_i-(l_cen0-l_p0),ldist_i] = z_sigmas0[np.where(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0)>max(np.mean(tpas_nc0[stime:etime,:,len_i,cent0+w_p0+l_dist0[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist0[ldist_i]]+w_col

for ldist_i in range(len(l_dist2)):
    for len_i in range(l_cen2-l_p2,l_cen2+l_p2):
        maxhe1[len_i-(l_cen2-l_p2),ldist_i] = z_sigmas1[np.where(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0)>max(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist2[ldist_i]] - z_sigmas1[np.where(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0)>max(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0))*perc)[0][0],len_i,l_dist2[ldist_i]]
        maxze1[len_i-(l_cen2-l_p2),ldist_i] = z_sigmas1[np.where(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0)>max(np.mean(tpas_nc1[stime:etime,:,len_i,cent2+w_p2+l_dist2[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist2[ldist_i]]+w_col

print('he1')
for ldist_i in range(len(l_dist1)):
    for len_i in range(l_cen1-l_p1,l_cen1+l_p1):
        maxhe2[len_i-(l_cen1-l_p1),ldist_i] = z_sigmas2[np.where(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0)>max(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist1[ldist_i]] - z_sigmas2[np.where(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0)>max(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0))*perc)[0][0],len_i,l_dist1[ldist_i]]
        maxhe3[len_i-(l_cen1-l_p1),ldist_i] = z_sigmas3[np.where(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0)>max(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist1[ldist_i]] - z_sigmas3[np.where(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0)>max(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0))*perc)[0][0],len_i,l_dist1[ldist_i]]
        maxze2[len_i-(l_cen1-l_p1),ldist_i] = z_sigmas2[np.where(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0)>max(np.mean(tpas_nc2[stime:etime,:,len_i,cent1+w_p1+l_dist1[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist1[ldist_i]]+w_col
        maxze3[len_i-(l_cen1-l_p1),ldist_i] = z_sigmas3[np.where(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0)>max(np.mean(tpas_nc3[stime:etime,:,len_i,cent1+w_p3+l_dist1[ldist_i]],axis=0))*perc)[0][-1],len_i,l_dist1[ldist_i]]+w_col


mean_ze = np.array([np.mean(maxze0,axis=0),np.mean(maxze1,axis=0),np.mean(maxze2,axis=0),np.mean(maxze3,axis=0)])/lb0
mean_he = np.array([np.mean(maxhe0,axis=0),np.mean(maxhe1,axis=0),np.mean(maxhe2,axis=0),np.mean(maxhe3,axis=0)])/lb0
'''

x_i_0 = l_real0/lb0
x_i_1 = l_real2/lb1
x_i_2 = l_real3/lb2
x_i_3 = l_real1/lb3
x_i_4 = l_real4/lb4
x_i_5 = l_real5/lb5
x_i_6 = l_real6/lb6

'''
############################
# vertical height at end of initial mixing region (dist_e = 500 m away from pipe)
############################
#ind_end0 = np.where(zm_plt0==max(zm_plt0))[0][0]
#ind_end1 = np.where(zm_plt1==max(zm_plt1))[0][0]
cx3 = 80
ind_end0 = -1
ind_end1 = -1
ind_end2 = -1
ind_end3 = (np.abs(x_i_3-cx3)).argmin()
ind_end4 = -1
# normalized
#zm_end = [zm_lb0[ind_end0],zm_lb1[ind_end1],zm_lb2[ind_end2],zm_lb3[ind_end3]]
# in m
#zm_end_real = np.array([zm0_arr[ind_end0],zm1_arr[ind_end1],zm2_arr[ind_end2],zm_lb3[ind_end3]]) 
#zm_end_real = np.array([zm0_arr[ind_end0],zm1_arr[ind_end1],zm2_arr[ind_end2],zm_lb3[ind_end3],zm_lb4[ind_end4]]) 
zm_end = [mean_zm[0][ind_end0],mean_zm[1][ind_end1],mean_zm[2][ind_end2],mean_zm[3][ind_end3]]
'''

############################################
# find length of inital mixing region
# by plotting dilution vs distance from pipe
############################################
# calculate minimum dilution (from max concentration)
# SqN/(b**(2./3))
S4 = np.empty((len(avg_prof4)))
S3 = np.empty((len(avg_prof3)))
S5 = np.empty((len(avg_prof5)))
S6 = np.empty((len(avg_prof6)))
for s_i in range(len(S4)):
    S4[s_i] = np.max(avg_prof4[s_i]) 

for s_i in range(len(S3)):
    S3[s_i] = np.max(avg_prof3[s_i]) 
    S5[s_i] = np.max(avg_prof5[s_i]) 
#    S4[s_i] = np.max(avg_prof4[s_i][1:]) # exclude first grid point of profile because
                                         # very high value there for F = 100

for s_i in range(len(S6)):
    S6[s_i] = np.max(avg_prof6[s_i]) 

'''
plt_S0 = ((tpas/S0)*(q*N0))/(buoy_roms**(2./3))
plt_S1 = ((tpas/S1)*(q*N1))/(buoy_roms**(2./3))
plt_S2 = ((tpas/S2)*(q*N2))/(buoy_roms**(2./3))
plt_S3 = ((tpas/S3)*(q*N3))/(buoy_roms**(2./3))
'''
plt_S4 = ((tpas/S4)*(q*N4))/(buoy_roms**(2./3))
plt_S2 = ((tpas/S3)*(q*N3))/(buoy_roms**(2./3))
plt_S5 = ((tpas/S5)*(q*N5))/(buoy_roms**(2./3))
plt_S6 = ((tpas/S6)*(q*N6))/(buoy_roms**(2./3))
'''

##########################################
# find thickness of established wastefield 
# (from previous calculation, dist_e = ~500 m away from pipe)
# at last time step
#########################################
dist_e0 = l_dist0[ind_end0] # distance away from pipe that represents established wastefield
dist_e1 = l_dist2[ind_end1]
dist_e2 = l_dist1[ind_end2]
dist_e3 = l_dist1[ind_end3]
#dist_e4 = l_dist1[ind_end4]
perc = 0.05
min_e0 = max(avg_prof0[ind_end0])*perc # minimum effluent detection is 5% of maximum conc
min_e1 = max(avg_prof1[ind_end1])*perc 
min_e2 = max(avg_prof2[ind_end2])*perc 
min_e3 = max(avg_prof3[ind_end3])*perc 
#min_e4 = max(avg_prof4[ind_end4])*perc 
# find where tpas averaged along pipe dist_e m away is greater than 5% of maximum conc
h0 = np.where(np.mean(np.mean(tpas_nc0[stime:etime,:,l_cen0-l_p0:l_cen0+l_p0,cent0+w_p0+dist_e0],axis=0),axis=1)>min_e0)[0]
h1 = np.where(np.mean(np.mean(tpas_nc1[stime:etime,:,l_cen2-l_p2:l_cen2+l_p2,cent2+w_p2+dist_e1],axis=0),axis=1)>min_e1)[0]
h2 = np.where(np.mean(np.mean(tpas_nc2[stime:etime,:,l_cen1-l_p1:l_cen1+l_p1,cent1+w_p1+dist_e2],axis=0),axis=1)>min_e2)[0]
h3 = np.where(np.mean(np.mean(tpas_nc3[stime:etime,:,l_cen1-l_p1:l_cen1+l_p1,cent1+w_p3+dist_e3],axis=0),axis=1)>min_e3)[0]
#h4 = np.where(np.mean(np.mean(tpas_nc4[stime:etime,:,l_cen1-l_p1:l_cen1+l_p1,cent1+w_p1+dist_e4],axis=0),axis=1)>min_e4)[0]
# find vertical depth at center of pipe dist_e m away
z_lev0 = z_sigmas0[h0,cent0,cent0+w_p0+dist_e0]
z_lev1 = z_sigmas1[h1,cent2,cent2+w_p2+dist_e1]
z_lev2 = z_sigmas2[h2,cent1,cent1+w_p1+dist_e2]
z_lev3 = z_sigmas3[h3,cent1,cent1+w_p3+dist_e3]
#z_lev4 = z_sigmas4[h4,cent1,cent1+w_p1+dist_e4]
# find thickness
he0 = (z_lev0[-1]-z_lev0[0])/lb0
he1 = (z_lev1[-1]-z_lev1[0])/lb1
he2 = (z_lev2[-1]-z_lev2[0])/lb2
he3 = (z_lev3[-1]-z_lev3[0])/lb3
#he4 = (z_lev4[-1]-z_lev4[0])/lb4

h_plt = np.array([he0,he1,he2,he3])

########################################################
# find minimum dilution at end of established wastefield
# (dist_e m away from pipe) averaged over time SmqN/b
########################################################
Sm_roms = np.array([plt_S0[ind_end0],plt_S1[ind_end1],plt_S2[ind_end2],plt_S3[ind_end3]])

###################################
# find top of established wastefield dist_e m away from pipe)
# from thickness calculation
###################################
ze = np.array([(z_lev0[-1]+w_col)/lb0,(z_lev1[-1]+w_col)/lb1,(z_lev2[-1]+w_col)/lb2,(z_lev3[-1]+w_col)/lb3])
'''


###########################
# RSB paper equations
# to plot next to my experiments
# and get data from rsb table
###############################
# Sm_rsb          RSB 1 Eq. 14
# Sm_rsb0         RSB 1 Eq. 12
# xi_lb_rsb       RSB 2 Eq. 11
# ze_lb_rsb_0     RSB 1 Eq. 16a
# ze_lb_rsb       RSB 1 Eq. 17
# zm_lb_rsb_0     RSB 1 Eq. 16c
# zm_lb_rsb       RSB 1 Eq. 18
# he_lb_rsb_0     RSB 1 Eq. 16b
# he_lb_rsb       RSb 1 Eq. 17

F_range_Sm = np.linspace(.1,100,1000)
F_range_Sm0 = np.linspace(0,.1,5)
F_range_1 = np.linspace(1,100,1000)
F_range_0 = np.linspace(0,.01,5) # between 0 and .1 F number 
F_range_discrete = np.array([0.1,1,10,100])


# equations

Sm_rsb = 2.19*(F_range_Sm**(1./6)) - 0.52 
Sm_rsb0 = np.full(5,0.97)

#xi_lb_rsb = 8.5*(F_range**(1./3))
xi_lb_rsb = 8.5*(F_range_discrete**(1./3)) # RSB2 eq 11


ze_lb_rsb_0 = np.full(5,2.6) 
ze_lb_rsb = 2.5*(F_range_1**(-1./6))


zm_lb_rsb_0 = np.full(5,1.7)
zm_lb_rsb = 1.5*(F_range_1**(-1./6))

he_lb_rsb_0 = np.full(5,1.8)
he_lb_rsb = 2.5*(F_range_1**(-1./6))

# rsb experiment parameters series 3 and 4 cross currents only
rsb_c = pd.read_csv('rsb_table_cross.csv',header=None)
#0                testnumber
#1    established_wastefield if 1, then parameters used for "end of initial mixing"/
#                            established wastefield
#2              port_spac_cm
#3      discharge_port_cm3_s cm^{3} s^{-1}
#4                 num_ports 50
#5                 g_prime_0 cm s^{-2}
#6                         N buoyancy frequency s^{-1}
#7                        Sa average dilution
#8                        Sm smallest dilution
#9    height_min_dilution_cm zm
#10  wastefield_thickness_cm he
#11 height_top_wastefield_cm ze
#12           sample_dist_cm xi
#13              current_dir
#14                   u_cm_s

# discharge per port
Qj = np.array(rsb_c[3][1:].astype(float))/(10**6) # convert to m3/s
# multiply by number of ports (50) and divide by length to get q
q_rsb_d = (Qj*float(rsb_c[4][1]))/1.2 # Q/L; m2/s
# buoyancy b = g'0*q
buoy_rsb_d = q_rsb_d*(np.array(rsb_c[5][1:].astype(float))/100) # m3/s3
u_rsb_d = np.array(rsb_c[14][1:].astype(float))/100
F_rsb_d = (u_rsb_d**3)/buoy_rsb_d
N_rsb_d = np.array(rsb_c[6][1:].astype(float))
lb_rsb_d = (buoy_rsb_d**(1./3))/N_rsb_d
# height and thickness converted to m and then nondimensionalized
zm_rsb_d = (np.array(rsb_c[9][1:].astype(float))/100.)/lb_rsb_d
he_rsb_d = (np.array(rsb_c[10][1:].astype(float))/100.)/lb_rsb_d
ze_rsb_d = (np.array(rsb_c[11][1:].astype(float))/100.)/lb_rsb_d
xi_d     = (np.array(rsb_c[12][1:].astype(float))/100.)/lb_rsb_d
# smallest dilutions
Sm_d = np.array(rsb_c[8][1:].astype(float))
Sm_d[Sm_d==0] = np.nan
SmqN_d = (Sm_d*q_rsb_d*N_rsb_d)/(buoy_rsb_d**(2./3))

# 0 values are no data in table
zm_rsb_d[zm_rsb_d==0] = np.nan
ze_rsb_d[ze_rsb_d==0] = np.nan
he_rsb_d[he_rsb_d==0] = np.nan

# separate by Froude numbers 0, .1, 1, 10 indices
F_rsb_0_d   = np.where(F_rsb_d==0)[0]
F_rsb_01_d  = np.where((F_rsb_d<.2)&(F_rsb_d>.05))[0]
F_rsb_1_d   = np.where((F_rsb_d<2)&(F_rsb_d>.5))[0]
F_rsb_10_d  = np.where((F_rsb_d<15)&(F_rsb_d>5))[0]
F_rsb_100_d = np.where(F_rsb_d>50)[0]

# separate by established wastefield for F vs SmqN, zm, ze, he calculations
est_w = np.where(np.array(rsb_c[1][1:].astype(float))==1)[0]
est0  = np.array(list(sorted(set(est_w) & set(F_rsb_0_d))))
est01 = np.array(list(sorted(set(est_w) & set(F_rsb_01_d))))
est1  = np.array(list(sorted(set(est_w) & set(F_rsb_1_d))))
est10 = np.array(list(sorted(set(est_w) & set(F_rsb_10_d))))
est100= np.array(list(sorted(set(est_w) & set(F_rsb_100_d))))

# eq. 11 RSB 2 closest approximation
Sm_11 = np.array([1,SmqN_d[F_rsb_1_d][1],SmqN_d[F_rsb_10_d][3],SmqN_d[F_rsb_100_d][1]])

'''
# same established wastefield curve for ROMS
#xi_lb_rom = [ind_end0,ind_end1,ind_end2,ind_end3,ind_end4]
#xi_lb_rom = np.load('F100_xi_lb_rom.npy')
Sm_xi_rom = np.load('Sm_xi_rom.npy')
# edit F = 100 dilution to 1 m result 
xi_lb_rom[4] = 53
Sm_xi_rom[4] = plt_S0[ind_end0]
np.save('F100_1m_x_end.npy',x_i_0[-1])
'''


plt_S0 = np.load('F100_1m_S.npy')
plt_S1 = np.load('F100_10m_S.npy')
plt_S3 = np.load('F100_5x300_S.npy')

######################
# plotting
#####################
#plt.ion()
axisfont = 16
tick_size=14

#lab0 = 'F = %.4f'%F_plt[0]+' 1 m'
#lab1 = 'F = %.4f'%F_plt[1]+' 10 m'
#lab2 = 'F = %.4f'%F_plt[2]
#lab3 = 'F = %.4f'%F_plt[3]
lab0 = 'F = 100 1 m'
lab1 = 'F = 100 10 m'
lab2 = 'F = 100 10 m KPP'
lab3 = 'F = 100 3 m 5x300'
lab4 = 'F = 100 10 m gray mix'
lab5 = 'F = 100 10 m no slip'
lab6 = 'F = 100 3 m gray mix'
rsb_lab = 'RSB'
rom_lab = 'ROMS'
lab_F0   = 'F = 1 m'
lab_F01  = 'F = 10 m'
lab_F1   = 'F = 3 m 10x300'
lab_F10  = 'F = 3 m 5x300'
lab_F100 = 'F = 100'
rsb_col_F0   = 'gray'
rsb_col_F01  = 'brown'
rsb_col_F1   = 'blue'
rsb_col_F10  = 'yellow'
rsb_col_F100 = 'limegreen'
rsb_col = 'red'
rsb_col0 = 'red'
rom_col = 'blue'
Fm0   = 'o'
Fm01  = '^'
Fm1   = 's'
Fm10  = 'x'
Fm100 = 'D'
F100_mark = 'P'

ltx = .01
lsty='-'
lw = 2

savename = 'exp_F100_vert_'

'''
# vertical height of min dilution vs distance from pipe
# plot zm/lb vs distance from diffuser
fig0,ax = plt.subplots(1,1,figsize=[12,8])
#ax.set_xlim([0,3])
#ax.set_ylim([0,2])
ax.plot(x_i_0,mean_zm[0],label=lab0,c=rsb_col_F0,linewidth=lw)
ax.plot(x_i_1,mean_zm[1],label=lab1,c=rsb_col_F01,linewidth=lw)
ax.plot(x_i_2,mean_zm[2],label=lab2,c=rsb_col_F1,linewidth=lw)
ax.plot(x_i_3,mean_zm[3],label=lab3,c=rsb_col_F10,linewidth=lw)
#ax.plot(x_i_4,zm_plt4,label=lab4,c=rsb_col_F100)
#ax.scatter(xi_d[F_rsb_0_d],zm_rsb_d[F_rsb_0_d],label=lab_F0,marker=Fm0,c=rsb_col_F0)
#ax.scatter(xi_d[F_rsb_01_d],zm_rsb_d[F_rsb_01_d],label=lab_F01,marker=Fm01,c=rsb_col_F01)
#ax.scatter(xi_d[F_rsb_1_d],zm_rsb_d[F_rsb_1_d],label=lab_F1,marker=Fm1,c=rsb_col_F1)
#ax.scatter(xi_d[F_rsb_10_d],zm_rsb_d[F_rsb_10_d],label=lab_F10,marker=Fm10,c=rsb_col_F10)
ax.scatter(xi_d[F_rsb_100_d],zm_rsb_d[F_rsb_100_d],label=lab_F100,marker=Fm100,c=rsb_col_F100)
ax.set_ylabel('z$_{\mathrm{m}}$/l$_{\mathrm{b}}$',fontsize=axisfont)
#ax.set_xlabel('Froude number, F = u$^3$/b',fontsize=axisfont)
ax.set_xlabel('distance from pipe in x-direction, x/l$_{\mathrm{b}}$',fontsize=axisfont)
ax.legend(loc='best')
ax.grid(True)
plt.savefig(savename+'z_vs_xdist.pdf',bbox_inches='tight')

# vertical height at end of initial mixing region dist_e m away from pipe)
fig1,ax = plt.subplots(1,1,figsize=[12,8])
ax.plot(F_plt,zm_end,label=rom_lab,linewidth=2)
ax.grid(True)
ax.scatter(100,zm_end[0],label=lab0,c=rsb_col_F0,marker=F100_mark)
ax.scatter(100,zm_end[1],label=lab1,c=rsb_col_F01,marker=F100_mark)
ax.scatter(100,zm_end[2],label=lab2,c=rsb_col_F1,marker=F100_mark)
ax.scatter(100,zm_end[3],label=lab3,c=rsb_col_F10,marker=F100_mark)
ax.plot(F_range_1,zm_lb_rsb,label=rsb_lab,color=rsb_col)
ax.plot(F_range_0,zm_lb_rsb_0,color=rsb_col)
#ax.scatter(F_rsb_d[est0],zm_rsb_d[est0],marker=Fm0,c=rsb_col_F0,label=lab_F0)
#ax.scatter(F_rsb_d[est01],zm_rsb_d[est01],marker=Fm01,c=rsb_col_F01,label=lab_F01)    
#ax.scatter(F_rsb_d[est1],zm_rsb_d[est1],marker=Fm1,c=rsb_col_F1,label=lab_F1)
#ax.scatter(F_rsb_d[est10],zm_rsb_d[est10],marker=Fm10,c=rsb_col_F10,label=lab_F10)    
ax.scatter(F_rsb_d[est100],zm_rsb_d[est100],marker=Fm100,c=rsb_col_F100,label=lab_F100)
ax.set_xlabel('Froude number, F = u$^3$/b',fontsize=axisfont)
ax.set_ylabel('z$_{\mathrm{m}}$/l$_{\mathrm{b}}$',fontsize=axisfont)
ax.set_xscale('symlog',basex=10,linthreshx=ltx)
ax.set_xlim(left=0)
xax = fig1.gca().xaxis
xax.set_minor_locator(MinorSymLogLocator(ltx))
ax.legend(loc='best')
plt.savefig(savename+'F_vs_zm.pdf',bbox_inches='tight') 
'''

# dilution vs distance from pipe to determine initial mixing region
fig2,ax = plt.subplots(1,1,figsize=[12,8])
ax.plot(x_i_0,plt_S0,label=lab0,c=rsb_col_F0,linewidth=lw) 
ax.plot(x_i_1,plt_S1,label=lab1,c=rsb_col_F01,linewidth=lw)
ax.plot(x_i_2,plt_S2,label=lab2,c=rsb_col_F1,linewidth=lw)
ax.plot(x_i_3,plt_S3,label=lab3,c=rsb_col_F10,linewidth=lw) 
ax.plot(x_i_4,plt_S4,label=lab4,c=rsb_col_F100,linewidth=lw) 
ax.plot(x_i_5,plt_S5,label=lab5,c='pink',linewidth=lw) 
ax.plot(x_i_6,plt_S6,label=lab6,c='purple',linewidth=lw) 
#ax.plot(xi_lb_rsb,Sm_11,c=rsb_col,label=rsb_lab,linestyle='--',linewidth=lw)
#ax.plot(xi_lb_rom,Sm_xi_rom,c='navy',label='ROMS',linestyle='--',linewidth=lw)
# don't plot first 4 values because plume doesn't start rising until 
# farther away from pipe
#ax.plot(x_i_4[167:],plt_S4[167:],label=lab4,c=rsb_col_F100) 
#ax.plot(x_i_4,plt_S4,label=lab4,c=rsb_col_F10)
#ax.scatter(xi_d[F_rsb_0_d],SmqN_d[F_rsb_0_d],label=lab_F0,marker=Fm0,c=rsb_col_F0)
#ax.scatter(xi_d[F_rsb_01_d],SmqN_d[F_rsb_01_d],label=lab_F01,marker=Fm01,c=rsb_col_F01)
#ax.scatter(xi_d[F_rsb_1_d],SmqN_d[F_rsb_1_d],label=lab_F1,marker=Fm1,c=rsb_col_F1)
#ax.scatter(xi_d[F_rsb_10_d],SmqN_d[F_rsb_10_d],label=lab_F10,marker=Fm10,c=rsb_col_F10)
ax.scatter(xi_d[F_rsb_100_d],SmqN_d[F_rsb_100_d],label=lab_F100,marker=Fm100,c=rsb_col_F100)
#ax.plot(xi_lb_rsb,label=rsb_lab)
ax.legend(loc='best')
ax.set_xlabel('distance from pipe in x-direction, x/l$_{\mathrm{b}}$',fontsize=axisfont)
ax.set_ylabel('SqN/b'+r'$^{\frac{2}{3}}$',fontsize=axisfont)
ax.tick_params(axis='both',which='major',labelsize=tick_size)
ax.grid(True)
plt.savefig(savename+'S_vs_xdist.png',bbox_inches='tight')
