#########################
# plot Wimp
# from random output
#########################

import os
import sys
sys.path.append('/data/project3/minnaho/global')
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import ROMS_depths as depths

plt.ion()

########
# 10 m
########

#print('pipes_10m_hydro_nogray_closed')
#datanc = Dataset('pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_dia.20000101054500.nc','r')
#gridnc = Dataset('pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_grd.nc','r')
#randnc = Dataset('pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_rnd.20000101054500.nc','r')
#
#cent = 65
#dx = 10.

#############
# test case
#############
datanc = Dataset('pipes_dt15_dx300_omega_Wimp_checkdiag_his.20000101000000.nc','r')
diagnc = Dataset('pipes_dt15_dx300_omega_Wimp_checkdiag_dia.20000101050000.nc','r')
gridnc = Dataset('pipes_dt15_dx300_omega_Wimp_checkdiag_grd.nc','r')
randnc = Dataset('pipes_dt15_dx300_omega_Wimp_checkdiag_rnd.20000101050000.nc','r')

dx = 300.
cent = 50

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

x_rho = np.array(gridnc.variables['x_rho']) # x locations of rho points
y_rho = np.array(gridnc.variables['y_rho']) # y locations of rho points

wimp = np.array(randnc.variables['Wimp'])/dx**2

figw = 14
figh = 8
axsize = 16
fig,ax = plt.subplots(1,1,figsize=[figw,figh],constrained_layout=True)
#stressplot = ax.pcolor(x_rho[cent,:],zw[:,cent,:],wimp[-1,:,cent,:],cmap='bwr',vmin=-0.2,vmax=0.2)
stressplot = ax.pcolor(x_rho[cent,:],zw[:,cent,:],wimp[-1,:,cent,:],cmap='bwr',vmin=-1E-18,vmax=1E-18)
cb0 = fig.colorbar(stressplot)
cb0.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb0.set_label(r'$\omega$ implicit m/s',fontsize=axsize)
ax.tick_params(axis='both',which='major',labelsize=axsize)
ax.set_ylabel('z (m)',fontsize=axsize)
ax.set_xlabel('x (m)',fontsize=axsize)




