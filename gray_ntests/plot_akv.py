import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import ROMS_depths as depths
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
plt.ion()

Zob = 1E-2
vonKar = 0.4

gridnc = Dataset('pipes_10m_hydro_gray_impnoslip_n3_grd.nc','r')

n0datanc = Dataset('pipes_10m_hydro_gray_impnoslip_n0_his.20000101000000.nc','r')
n1datanc = Dataset('pipes_10m_hydro_gray_impnoslip_n1_his.20000101000000.nc','r')
n2datanc = Dataset('pipes_10m_hydro_gray_impnoslip_n2_his.20000101000000.nc','r')
n3datanc = Dataset('pipes_10m_hydro_gray_impnoslip_n3_his.20000101000000.nc','r')
n4datanc = Dataset('pipes_10m_hydro_gray_impnoslip_n4_his.20000101000000.nc','r')
n5datanc = Dataset('pipes_10m_hydro_gray_impnoslip_n5_his.20000101000000.nc','r')
RIdatanc = Dataset('pipes_10m_hydro_RImix_noslip_his.20000101000000.nc','r')
n2newdatanc = Dataset('pipes_10m_hydro_gray_impnoslip_n2_newhn_his.20000101000000.nc','r')

un0 = np.array(n0datanc.variables['u'])
un1 = np.array(n1datanc.variables['u'])
un2 = np.array(n2datanc.variables['u'])
un3 = np.array(n3datanc.variables['u'])
un4 = np.array(n4datanc.variables['u'])
un5 = np.array(n5datanc.variables['u'])
uRI = np.array(RIdatanc.variables['u'])
un2new = np.array(n2newdatanc.variables['u'])

hnc = np.array(gridnc.variables['h'])

zwn0 = depths.get_zw_zeta_tind(n0datanc,gridnc,-1)
zwn1 = depths.get_zw_zeta_tind(n1datanc,gridnc,-1)
zwn2 = depths.get_zw_zeta_tind(n2datanc,gridnc,-1)
zwn3 = depths.get_zw_zeta_tind(n3datanc,gridnc,-1)
zwn4 = depths.get_zw_zeta_tind(n4datanc,gridnc,-1)
zwn5 = depths.get_zw_zeta_tind(n5datanc,gridnc,-1)
zwRI = depths.get_zw_zeta_tind(RIdatanc,gridnc,-1)
zwn2new = depths.get_zw_zeta_tind(n2newdatanc,gridnc,-1)

zrn0 = depths.get_zr_zeta_tind(n0datanc,gridnc,-1)
zrn1 = depths.get_zr_zeta_tind(n1datanc,gridnc,-1)
zrn2 = depths.get_zr_zeta_tind(n2datanc,gridnc,-1)
zrn3 = depths.get_zr_zeta_tind(n3datanc,gridnc,-1)
zrn4 = depths.get_zr_zeta_tind(n4datanc,gridnc,-1)
zrn5 = depths.get_zr_zeta_tind(n5datanc,gridnc,-1)
zrRI = depths.get_zr_zeta_tind(RIdatanc,gridnc,-1)
zrn2new = depths.get_zr_zeta_tind(n2newdatanc,gridnc,-1)

akvn0 = np.array(n0datanc.variables['Akv'])[-1]
akvn1 = np.array(n1datanc.variables['Akv'])[-1]
akvn2 = np.array(n2datanc.variables['Akv'])[-1]
akvn3 = np.array(n3datanc.variables['Akv'])[-1]
akvn4 = np.array(n4datanc.variables['Akv'])[-1]
akvn5 = np.array(n5datanc.variables['Akv'])[-1]
akvRI = np.array(RIdatanc.variables['Akv'])[-1]
akvn2new = np.array(n2newdatanc.variables['Akv'])[-1]


fig,ax = plt.subplots(1,1)
ax.plot(akvn0[:,65,90],zwn0[:,65,90],label='n = 0')
ax.plot(akvn1[:,65,90],zwn1[:,65,90],label='n = 1')
ax.plot(akvn2[:,65,90],zwn2[:,65,90],label='n = 2')
ax.plot(akvn3[:,65,90],zwn3[:,65,90],label='n = 3')
ax.plot(akvn4[:,65,90],zwn4[:,65,90],label='n = 4')
ax.plot(akvn5[:,65,90],zwn5[:,65,90],label='n = 5')
ax.plot(akvRI[:,65,90],zwRI[:,65,90],label='RImix')
ax.plot(akvn2new[:,65,90],zwn2new[:,65,90],label='n = 2')

ax.scatter(akvn0[:,65,90],zwn0[:,65,90],s=8)
ax.scatter(akvn1[:,65,90],zwn1[:,65,90],s=8)
ax.scatter(akvn2[:,65,90],zwn2[:,65,90],s=8)
ax.scatter(akvn3[:,65,90],zwn3[:,65,90],s=8)
ax.scatter(akvn4[:,65,90],zwn4[:,65,90],s=8)
ax.scatter(akvn5[:,65,90],zwn5[:,65,90],s=8)
ax.scatter(akvRI[:,65,90],zwRI[:,65,90],s=8)
ax.scatter(akvn2new[:,65,90],zwn2new[:,65,90],s=8)

#ax.set_ylim([-60,-30])
ax.set_xlabel('Akv')
ax.set_ylabel('depth')
ax.legend()

fig1,ax1 = plt.subplots(1,1)
ax1.plot(un0[-1,:,65,90],zrn0[:,65,90],label='n = 0')
ax1.plot(un1[-1,:,65,90],zrn1[:,65,90],label='n = 1')
ax1.plot(un2[-1,:,65,90],zrn2[:,65,90],linestyle='--',label='n = 2')
ax1.plot(un3[-1,:,65,90],zrn3[:,65,90],label='n = 3')
ax1.plot(un4[-1,:,65,90],zrn4[:,65,90],label='n = 4')
ax1.plot(un5[-1,:,65,90],zrn5[:,65,90],label='n = 5')
ax1.plot(uRI[-1,:,65,90],zrRI[:,65,90],linestyle='--',label='RImix')
ax1.plot(un2new[-1,:,65,90],zrn2new[:,65,90],label='n = 2 new')

ax1.scatter(un0[-1,:,65,90],zrn0[:,65,90],s=4)
ax1.scatter(un1[-1,:,65,90],zrn1[:,65,90],s=4)
ax1.scatter(un2[-1,:,65,90],zrn2[:,65,90],s=4)
ax1.scatter(un3[-1,:,65,90],zrn3[:,65,90],s=4)
ax1.scatter(un4[-1,:,65,90],zrn4[:,65,90],s=4)
ax1.scatter(un5[-1,:,65,90],zrn5[:,65,90],s=4)
ax1.scatter(uRI[-1,:,65,90],zrRI[:,65,90],s=4)
ax1.scatter(un2new[-1,:,65,90],zrn2new[:,65,90],s=4)

ax1.set_xlabel('u (m/s)')
ax1.set_ylabel('depth')
ax1.legend()

fig2,ax2 = plt.subplots(1,1)
ax2.plot(un0[-1,:,65,90],zrn0[:,65,90],label='n = 0')
ax2.plot(un1[-1,:,65,90],zrn1[:,65,90],label='n = 1')
ax2.plot(un2[-1,:,65,90],zrn2[:,65,90],linestyle='--',label='n = 2')
ax2.plot(un3[-1,:,65,90],zrn3[:,65,90],label='n = 3')
ax2.plot(un4[-1,:,65,90],zrn4[:,65,90],label='n = 4')
ax2.plot(un5[-1,:,65,90],zrn5[:,65,90],label='n = 5')
ax2.plot(uRI[-1,:,65,90],zrRI[:,65,90],linestyle='--',label='RImix')
ax2.plot(un2new[-1,:,65,90],zrn2new[:,65,90],label='n = 2 new')

ax2.scatter(un0[-1,:,65,90],zrn0[:,65,90],s=4)
ax2.scatter(un1[-1,:,65,90],zrn1[:,65,90],s=4)
ax2.scatter(un2[-1,:,65,90],zrn2[:,65,90],s=4)
ax2.scatter(un3[-1,:,65,90],zrn3[:,65,90],s=4)
ax2.scatter(un4[-1,:,65,90],zrn4[:,65,90],s=4)
ax2.scatter(un5[-1,:,65,90],zrn5[:,65,90],s=4)
ax2.scatter(uRI[-1,:,65,90],zrRI[:,65,90],s=4)
ax2.scatter(un2new[-1,:,65,90],zrn2new[:,65,90],s=4)

ax2.set_ylim([-60,-45])
ax2.set_xlim([0.02,0.065])
ax2.set_xlabel('u (m/s)')
ax2.set_ylabel('depth')
ax2.legend()
