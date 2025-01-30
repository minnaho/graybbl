#########################
# Calculate bottom stress

# tau = rho C_d u_b^2
# tau = stress
# rho = density
# C_d = drag coefficient
# u_b = u at the bottom

# C_d = k^2/(ln(zb/z0))^2
# zb = z near bottom (height of grid cell at bottom)
# z0 = roughness length (Zob)
#########################

import os
import sys
sys.path.append('/data/project3/minnaho/global')
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import ROMS_depths as depths

plt.ion()


vonKar = 0.41
Zob = 1E-2

########
# 10 m
########

print('pipes_10m_nh_grayn1')
datanc = Dataset('pipes_10m_nh_grayn1_his.20000101000000.nc','r')
diagnc = Dataset('pipes_10m_nh_grayn1_dia.20000101001500.nc','r')
gridnc = Dataset('pipes_10m_nh_grayn1_grd.nc','r')

#print('pipes_10m_nh_grayn2')
#datanc = Dataset('pipes_10m_nh_grayn2_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_10m_nh_grayn2_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_10m_nh_grayn2_grd.nc','r')

########
# 3 m 
########

#print('pipes_3m_dt1_nh_grayn1')
#datanc = Dataset('pipes_3m_dt1_nh_grayn1_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_3m_dt1_nh_grayn1_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_3m_dt1_nh_grayn1_grd.nc','r')

#print('pipes_3m_dt1_nh_grayn2')
#datanc = Dataset('pipes_3m_dt1_nh_grayn2_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_3m_dt1_nh_grayn2_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_3m_dt1_nh_grayn2_grd.nc','r')
#gridnc10 = Dataset('pipes_10m_nh_gray_n1_diag_grd.nc','r')

zw = depths.get_zw_zeta_tind(datanc,gridnc,-1)

# get height of zw's
Hz = np.diff(zw,axis=0)

C_d = (vonKar/np.log(1+0.5*Hz[0,:,:]/Zob))**2
u_b = np.array(datanc.variables['u'])[-1,0,:,:] # bottom u
v_b = np.array(datanc.variables['v'])[-1,0,:,:] # bottom v
rho = np.array(datanc.variables['rho'])[-1,0,:,:] # bottom density
x_rho = np.array(gridnc.variables['x_rho']) # x locations of rho points
y_rho = np.array(gridnc.variables['y_rho']) # y locations of rho points

u_b_rho = np.empty((u_b.shape[0],u_b.shape[1]-1))*np.nan
v_b_rho = np.empty((v_b.shape[0]-1,v_b.shape[1]))*np.nan

# move u_b to rho points
for i in range(u_b.shape[1]-1):
    u_b_rho[:,i] = 0.5*(u_b[:,i]+u_b[:,i+1])

tau = rho[:,1:-1]*C_d[:,1:-1]*u_b_rho

# find where x_rho is the same as 10 m resolution
# only for 3 m resolution!

#x_rho10 = np.array(gridnc10.variables['x_rho']) # x locations of rho points

#x_where = np.where(x_rho[0]<x_rho10[0])

figw = 14
figh = 8
axsize = 16
fig,ax = plt.subplots(1,1,figsize=[figw,figh],constrained_layout=True)
stressplot = ax.pcolor(x_rho[:,1:-1],y_rho[:,1:-1],tau,cmap='bwr',vmin=-0.006,vmax=0.006)
cb0 = fig.colorbar(stressplot)
cb0.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb0.set_label(r'$\tau_b$'+' (kg m$^{-1}$ s$^{-2}$)',fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)
ax.set_ylabel('y (m)',fontsize=axsize)
ax.set_xlabel('x (m)',fontsize=axsize)




