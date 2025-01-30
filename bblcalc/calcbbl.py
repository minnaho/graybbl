##############
# Calculate bottom boundary layer (BBL) thickness
# based on theory from 
# Boundary Layer Theory 9th edition by Schlichting & Gersten 
# compare ROMS gray mix BBL thickness vs theory
###############

import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
from netCDF4 import Dataset
import ROMS_depths as depths

vonKar = 0.41
U = 0.67 # velocity away from the boundary

nu = U/12 # probably wrong, from UP3 hyperviscosity

data10 = Dataset('/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_gray_max_his.20000101000000.nc','r')
grd10 = Dataset('/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_gray_max_grd.nc','r')

data3  = Dataset('/data/project9/minnaho/graybbl/F100_tests/pipe_3m_dt1_nh_F100_gray_nopipe_his.20000101000000.nc','r')
grd3  = Dataset('/data/project9/minnaho/graybbl/F100_tests/pipe_3m_dt1_nh_F100_gray_nopipe_grd.nc','r')

u10 = np.array(data10.variables['u'][-1,:,:,:])
u3  = np.array(data3.variables['u'][-1,:,:,:])

# halfway into x domain
x10 = 128
x3  = 512

y10 = 65
y3 = 257

# calculate nu (kinematic viscosity), nu = dynamic viscosity/density

# Reynolds number (with scaling using "plate length" l)

# Reynolds number x distance away from start of turbulent
# boundary layer
Rex10 = U*x10/nu
Rex3 = U*x3/nu


# theoretical turbulent BBL thickness delta
#theod = ((0.14*nu*Rex*np.log(Re))/(U*np.log(Rex)*vonKar))*((Cd/2)**0.5)
#theotd10 = 0.21*(x10/(np.log(U*x10/nu)))
#theotd3  = 0.21*(x3/(np.log(U*x3/nu)))
#
## theoretical laminar BBL thickness 
#theold10 = 5*(nu*x10/U)**0.5
#theold3 = 5*(nu*x3/U)**0.5

# theoretical according to 
# fluid-mechanics-fundamentals-and-applications-3rd-edition-
# cengel-and-cimbala-2014
theotd10 = 0.38*x10/(Rex10**0.2)
theotd3 = 0.38*x3/(Rex3**0.2)

z10 = depths.get_zw_zeta_tind(data10,grd10,-1)
z3  = depths.get_zw_zeta_tind(data3,grd3,-1)

dz10 = np.diff(z10[:,y10,x10])
dz3 = np.diff(z3[:,y3,x3])

# displacement distance delta1
smax10 = np.where(u10[:,y10,x10]==np.nanmax(u10[:,y10,x10]))[0][0]
smax3  = np.where(u3[:,y3,x3]==np.nanmax(u3[:,y3,x3]))[0][0]


delta1_10 = (1/U)*np.nansum((U-u10[:,y10,x10])*dz10)
delta1_3  = (1/U)*np.nansum((U-u3[:,y3,x3])*dz3)
