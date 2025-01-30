import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import ROMS_depths as depths
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

plt.ion()

data = Dataset('fb_computerd_graynew_minus1_his.20000101000000.nc','r')
grid = Dataset('fb_computerd_graynew_minus1_grd.nc','r')

zr = depths.get_zr_zeta_tind(data,grid,-1)
zw = depths.get_zw_zeta_tind(data,grid,-1)

akv = np.array(data.variables['Akv'][-1,:,0,0])
u = np.array(data.variables['u'][-1,:,0,0])

axsize = 16

fig1,ax1 = plt.subplots(1,1,figsize=[5,12])
ax1.plot(akv,zw[:,0,0])
ax1.scatter(akv,zw[:,0,0])
ax1.set_xlabel('Akv (m$^2$/s)',fontsize=axsize)
ax1.set_ylabel('depth (m)',fontsize=axsize)


fig2,ax2 = plt.subplots(1,1,figsize=[5,12])
ax2.plot(u,zr[:,0,0])
ax2.scatter(u,zr[:,0,0])
ax2.set_xlabel('u (m/s)',fontsize=axsize)
ax2.set_ylabel('depth (m)',fontsize=axsize)


ax1.tick_params(axis='both',which='major',labelsize=axsize)
ax2.tick_params(axis='both',which='major',labelsize=axsize)

