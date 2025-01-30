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

datanc = Dataset('pipes_10m_hydro_grayexchange1_his.20000101000000.0.nc','r')
datanc4 = Dataset('pipes_10m_hydro_gray4_his.20000101000000.0.nc','r')
gridnc = Dataset('pipes_10m_hydro_grayexchange1_grd.0.nc','r')
#nograync = Dataset('pipes_10m_hydro_1_his.20000101000000.0.nc','r')
nograync = Dataset('pipes_10m_hydro_freeslip_his.20000101000000.0.nc','r')

impbothnc = Dataset('pipes_10m_hydro_impboth_his.20000101000000.nc','r')
rimix_impbothnc = Dataset('pipes_10m_hydro_impboth_RImix_his.20000101000000.nc','r')
rimixnc = Dataset('pipes_10m_hydro_RImix_noimp_his.20000101000000.nc','r')
grayimpbothnc = Dataset('pipes_10m_hydro_impboth_gray_his.20000101000000.nc','r')

grayimpbothnc_4 = Dataset('pipes_10m_hydro_gray_impboth_his.20000101000000.0.nc','r')
grayimpnoslipnc = Dataset('pipes_10m_hydro_gray_impnoslip_his.20000101000000.0.nc','r')



uncnogray = np.array(nograync.variables['u'])
unc4 = np.array(datanc4.variables['u'])

unc_impboth = np.array(impbothnc.variables['u'])
unc_rimix_impboth = np.array(rimix_impbothnc.variables['u'])

unc_rimix = np.array(rimixnc.variables['u'])
unc_grayimpboth = np.array(grayimpbothnc.variables['u'])
unc_grayimpboth_4 = np.array(grayimpbothnc_4.variables['u'])
unc_grayimpnoslip = np.array(grayimpnoslipnc.variables['u'])


hnc = np.array(gridnc.variables['h'])

unc = np.array(datanc.variables['u'])
vnc = np.array(datanc.variables['v'])

cff = np.empty((130,258))

for i in range(unc.shape[3]-1):
    for j in range(vnc.shape[2]-1):
        cff[j,i] = np.sqrt(0.333333333333*(
                     (unc[-1,0,j,i]**2+unc[-1,0,j,i+1]**2)
                    +(unc[-1,0,j,i]*unc[-1,0,j,i+1])
                    +(vnc[-1,0,j,i]**2+vnc[-1,0,j+1,i]**2)
                    +(vnc[-1,0,j,i]*vnc[-1,0,j+1,i])
                                     ))

zw = depths.get_zw_zeta_tind(datanc,gridnc,-1)
zw4 = depths.get_zw_zeta_tind(datanc4,gridnc,-1)
zw_rimix = depths.get_zw_zeta_tind(rimixnc,gridnc,-1)
zw_rimix_impboth = depths.get_zw_zeta_tind(rimix_impbothnc,gridnc,-1)
zw_grayimpboth = depths.get_zw_zeta_tind(grayimpbothnc,gridnc,-1)
zw_grayimpboth_4 = depths.get_zw_zeta_tind(grayimpbothnc_4,gridnc,-1)
zw_grayimpnoslip = depths.get_zw_zeta_tind(grayimpnoslipnc,gridnc,-1)
zr = depths.get_zr_zeta_tind(datanc,gridnc,-1)
zr4 = depths.get_zr_zeta_tind(datanc4,gridnc,-1)
zr_freeslip = depths.get_zr_zeta_tind(nograync,gridnc,-1)
zr_impboth = depths.get_zr_zeta_tind(impbothnc,gridnc,-1)
zr_rimix_impboth = depths.get_zr_zeta_tind(rimix_impbothnc,gridnc,-1)
zr_rimix = depths.get_zr_zeta_tind(rimixnc,gridnc,-1)
zr_grayimpboth = depths.get_zr_zeta_tind(grayimpbothnc,gridnc,-1)
zr_grayimpboth_4 = depths.get_zr_zeta_tind(grayimpbothnc_4,gridnc,-1)
zr_grayimpnoslip = depths.get_zr_zeta_tind(grayimpnoslipnc,gridnc,-1)
# get height of zw's
Hz = np.diff(zw,axis=0)

rd = cff*(vonKar/np.log(1+0.5*Hz[0,:,:]/Zob))**2

hn = zw[5,:,:]+hnc

ustar = np.empty((130,258))
for i in range(unc.shape[3]-1):
    for j in range(vnc.shape[2]-1):
        ustar[j,i]=np.sqrt( 
                   np.sqrt((rd[j,i]*(unc[-1,0,j,i]*unc[-1,0,j,i+1])/2)**2+
                           (rd[j,i]*(vnc[-1,0,j,i]*vnc[-1,0,j+1,i])/2)**2))

ustar = ustar*40

akv_ana = ustar*(zw+hnc)/(vonKar*(1+((zw+hnc)/hn))**2)
akv_ana4 = ustar*(zw+hnc)/(vonKar*(1+((zw+hnc)/hn))**4)



akv_nc = np.array(datanc.variables['Akv'])[-1]
akv_nc4 = np.array(datanc4.variables['Akv'])[-1]
akv_nc_rimix = np.array(rimixnc.variables['Akv'])[-1]
akv_nc_rimix_impboth = np.array(rimix_impbothnc.variables['Akv'])[-1]
akv_nc_gray_impboth = np.array(grayimpbothnc.variables['Akv'])[-1]
akv_nc_gray_impboth_4 = np.array(grayimpbothnc_4.variables['Akv'])[-1]
akv_nc_gray_impnoslip = np.array(grayimpnoslipnc.variables['Akv'])[-1]

fig,ax = plt.subplots(1,1)
ax.plot(akv_nc[:,65,90],zw[:,65,90],label='roms gray power 2')
ax.plot(akv_nc4[:,65,90],zw4[:,65,90],label='roms gray power 4')
ax.plot(akv_nc_rimix[:,65,90],zw_rimix[:,65,90],label='roms RImix')
ax.plot(akv_nc_rimix_impboth[:,65,90],zw_rimix_impboth[:,65,90],label='roms RImix impboth')
ax.plot(akv_nc_gray_impboth[:,65,90],zw_grayimpboth[:,65,90],label='roms gray power 2 impboth')
ax.plot(akv_nc_gray_impboth_4[:,65,90],zw_grayimpboth_4[:,65,90],label='roms gray power 4 impboth')
ax.plot(akv_nc_gray_impnoslip[:,65,90],zw_grayimpnoslip[:,65,90],label='roms gray power 4 impboth')
ax.plot(akv_ana[:,65,90],zw[:,65,90],label='analytical power 2')
ax.plot(akv_ana4[:,65,90],zw[:,65,90],label='analytical power 4')
#ax.plot(akv_ana[:,65,90]*40,zw[:,65,90],label='analytical')
#ax.plot(akv_ana4[:,65,90]*40,zw[:,65,90],label='analytical power 4')
ax.set_xlabel('Akv')
ax.set_ylabel('depth')
ax.legend()

fig1,ax1 = plt.subplots(1,1)
ax1.plot(unc[-1,:,65,90],zr[:,65,90],label='gray power 2')
ax1.plot(unc4[-1,:,65,90],zr4[:,65,90],label='gray power 4')
ax1.plot(uncnogray[-1,:,65,90],zr_freeslip[:,65,90],label='free slip',linestyle='--')
#ax1.plot(unc_impboth[-1,:,65,90],zr_impboth[:,65,90],label='imp both')
#ax1.scatter(unc_impboth[-1,:,65,90],zr_impboth[:,65,90],label='imp both')
ax1.plot(unc_rimix[-1,:,65,90],zr_rimix[:,65,90],label='RImix',linestyle='-.')
ax1.plot(unc_rimix_impboth[-1,:,65,90],zr_rimix_impboth[:,65,90],label='RImix imp both',linestyle=':')
ax1.plot(unc_grayimpboth[-1,:,65,90],zr_grayimpboth[:,65,90],label='gray power 2 imp both')
ax1.plot(unc_grayimpboth_4[-1,:,65,90],zr_grayimpboth_4[:,65,90],label='gray power 4 imp both')
ax1.plot(unc_grayimpnoslip[-1,:,65,90],zr_grayimpnoslip[:,65,90],label='gray power 4 imp noslip')
ax1.set_xlabel('u (m/s)')
ax1.set_ylabel('depth')
ax1.legend()
