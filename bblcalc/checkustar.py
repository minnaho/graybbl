import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
import ROMS_depths as depths
from netCDF4 import Dataset

u0 = 0.67
v0 = 0
vonKar = 0.41
Zob = 1.0E-2

cff0=np.sqrt( 0.333333333333*(
      u0**2 +u0**2
             +u0*u0
       +v0**2+v0**2
             +v0*v0))


his = Dataset('fb_computerd_graynew_minus1_his.20000101000000.nc','r')
grd = Dataset('fb_computerd_graynew_minus1_grd.nc','r')
rnd = Dataset('fb_computerd_graynew_minus1_rnd.20000101000004.nc','r')

zw0 = depths.get_zw_zeta_tind(his,grd,0)

Hz0 = zw0[1,0,0]-zw0[0,0,0]

#rd0= cff0*(vonKar/np.log(1+0.5*Hz0/Zob))**2
rd0 = np.array(rnd.variables['rd_out'][0,0,0])

ustar0=np.sqrt(
      np.sqrt((rd0*(u0+u0)/2)**2+
              (rd0*(v0+v0)/2)**2))

u1 = np.array(his.variables['u'][-1,0,0,0])
v1 = np.array(his.variables['v'][-1,0,0,0])

cff1=np.sqrt( 0.333333333333*(
      u1**2 +u1**2
             +u1*u1
       +v1**2+v1**2
             +v1*v1))

zw1 = depths.get_zw_zeta_tind(his,grd,1)

Hz1 = zw1[1,0,0]-zw1[0,0,0]

#rd1= cff1*(vonKar/np.log(1+0.5*Hz1/Zob))**2
rd1 = np.array(rnd.variables['rd_out'][1,0,0])

ustar1=np.sqrt(
      np.sqrt((rd1*(u1+u1)/2)**2+
              (rd1*(v1+v1)/2)**2))

