import os
import sys
sys.path.append('/data/project3/minnaho/global')
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import ROMS_depths as depths
plt.ion()

# F = 100 (67 cm/s) case 10 m resolution (but diag vmix is wrong)

#datanc = Dataset('pipes_10m_nh_gray_impnoslip_n1_F100_diag_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_10m_nh_gray_impnoslip_n1_F100_diag_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_10m_nh_gray_impnoslip_n1_F100_diag_grd.nc','r')

###########
# test cases
###########
#cent = int(100/2)

#datanc = Dataset('pipes_dt60_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt60_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt60_grd.nc','r')

#datanc = Dataset('pipes_dt30_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt30_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt30_grd.nc','r')

#datanc = Dataset('pipes_dt15_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt15_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt15_grd.nc','r')

#datanc = Dataset('pipes_dt15_dx30_omega_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt15_dx30_omega_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt15_dx30_omega_grd.nc','r')

#datanc = Dataset('pipes_dt05_dx30_omega_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt05_dx30_omega_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt05_dx30_omega_grd.nc','r')

#datanc = Dataset('pipes_dt60_dx300_omega_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt60_dx300_omega_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt60_dx300_omega_grd.nc','r')

#datanc = Dataset('pipes_dt15_dx300_omega_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_dt15_dx300_omega_dia.20000101050000.nc','r')
#gridnc = Dataset('pipes_dt15_dx300_omega_grd.nc','r')

datanc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_dt15_dx300_omega_Wimp_checkdiag_his.20000101000000.nc','r')
diagnc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_dt15_dx300_omega_Wimp_checkdiag_dia.20000101050000.nc','r')
gridnc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_dt15_dx300_omega_Wimp_checkdiag_grd.nc','r')
randnc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_dt15_dx300_omega_Wimp_checkdiag_rnd.20000101050000.nc','r')

dx = 300.
cent = 50


###########
# 10 m 
###########

#cent = int(128/2)
#
#datanc = Dataset('pipes_10m_nh_grayn1_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_10m_nh_grayn1_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_10m_nh_grayn1_grd.nc','r')

#datanc = Dataset('pipes_10m_nh_grayn2_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_10m_nh_grayn2_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_10m_nh_grayn2_grd.nc','r')

#datanc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_his.20000101000000.nc','r')
#diagnc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_dia.20000101054500.nc','r')
#gridnc = Dataset('/data/project9/minnaho/graybbl/check_diag/pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_grd.nc','r')


###########
# 3 m 
###########

#cent = int(512/2)

#print('pipes_3m_dt1_nh_grayn1')
#datanc = Dataset('pipes_3m_dt1_nh_grayn1_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_3m_dt1_nh_grayn1_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_3m_dt1_nh_grayn1_grd.nc','r')

#print('pipes_3m_dt1_nh_grayn2')
#datanc = Dataset('pipes_3m_dt1_nh_grayn2_his.20000101000000.nc','r')
#diagnc = Dataset('pipes_3m_dt1_nh_grayn2_dia.20000101001500.nc','r')
#gridnc = Dataset('pipes_3m_dt1_nh_grayn2_grd.nc','r')

#################################################################
# begin reading and calculation
#################################################################

xslice_m = np.array(gridnc.variables['x_rho'])[cent]

udiss = np.array(diagnc.variables['u_dis'])[-1,:,:,:]
vdiss = np.array(diagnc.variables['v_dis'])[-1,:,:,:]

uvmix = np.array(diagnc.variables['u_vmx'])[-1,:,:,:]
vvmix = np.array(diagnc.variables['v_vmx'])[-1,:,:,:]

#upgr = np.array(diagnc.variables['u_pgr'])
#vpgr = np.array(diagnc.variables['v_pgr'])
#
#uadv = np.array(diagnc.variables['u_adv'])
#vadv = np.array(diagnc.variables['v_adv'])

unc = np.array(datanc.variables['u'])[-1,:,:,:]
vnc = np.array(datanc.variables['v'])[-1,:,:,:]

akv = np.array(datanc.variables['Akv'])[-1,:,:,:]

udiss[udiss>1E10] = np.nan
vdiss[vdiss>1E10] = np.nan

uvmix[uvmix>1E10] = np.nan
vvmix[vvmix>1E10] = np.nan

#upgr[upgr>1E10] = np.nan
#vpgr[vpgr>1E10] = np.nan
#
#uadv[uadv>1E10] = np.nan
#vadv[vadv>1E10] = np.nan

unc[unc>1E10] = np.nan
vnc[vnc>1E10] = np.nan


# depths and Hz at last time step
#zw = depths.get_zw_zeta_tind(datanc,gridnc,-1) # depths at w points
#zr = depths.get_zr_zeta_tind(datanc,gridnc,-1) # depths at rho points
zw = depths.get_zw_zeta_tind(datanc,gridnc,0) # depths at w points
zr = depths.get_zr_zeta_tind(datanc,gridnc,0) # depths at rho points
Hz = np.diff(zw,axis=0) # height of dz

# rho point variables to calculate
'''
udiss_rho = np.empty((udiss.shape[0],udiss.shape[1],vdiss.shape[2]-1,udiss.shape[3]-1))
vdiss_rho = np.empty((udiss.shape[0],udiss.shape[1],vdiss.shape[2]-1,udiss.shape[3]-1))

uvmix_rho = np.empty((uvmix.shape[0],uvmix.shape[1],vvmix.shape[2]-1,uvmix.shape[3]-1))
vvmix_rho = np.empty((uvmix.shape[0],uvmix.shape[1],vvmix.shape[2]-1,uvmix.shape[3]-1))

#upgr_rho = np.empty((upgr.shape[0],upgr.shape[1],vpgr.shape[2]-1,upgr.shape[3]-1))
#vpgr_rho = np.empty((upgr.shape[0],upgr.shape[1],vpgr.shape[2]-1,upgr.shape[3]-1))
#
#uadv_rho = np.empty((uadv.shape[0],uadv.shape[1],vadv.shape[2]-1,uadv.shape[3]-1))
#vadv_rho = np.empty((uadv.shape[0],uadv.shape[1],vadv.shape[2]-1,uadv.shape[3]-1))

u_rho = np.empty((unc.shape[0],unc.shape[1],vnc.shape[2]-1,unc.shape[3]-1))*np.nan
v_rho = np.empty((vnc.shape[0],vnc.shape[1],vnc.shape[2]-1,unc.shape[3]-1))*np.nan
'''

#for i in range(udiss.shape[3]-1):
#    print('calc '+str(i)+' of '+str(udiss.shape[3]-1))
#    for j in range(vdiss.shape[2]-1):
#        udiss_rho[:,:,j,i] = (udiss[:,:,j,i]+udiss[:,:,j,i+1])/2
#        vdiss_rho[:,:,j,i] = (vdiss[:,:,j,i]+vdiss[:,:,j+1,i])/2
#
#        uvmix_rho[:,:,j,i] = (uvmix[:,:,j,i]+uvmix[:,:,j,i+1])/2
#        vvmix_rho[:,:,j,i] = (vvmix[:,:,j,i]+vvmix[:,:,j+1,i])/2
#
#        upgr_rho[:,:,j,i] = (upgr[:,:,j,i]+upgr[:,:,j,i+1])/2
#        vpgr_rho[:,:,j,i] = (vpgr[:,:,j,i]+vpgr[:,:,j+1,i])/2
#
#        uadv_rho[:,:,j,i] = (uadv[:,:,j,i]+uadv[:,:,j,i+1])/2
#        vadv_rho[:,:,j,i] = (vadv[:,:,j,i]+vadv[:,:,j+1,i])/2
#
#        # calculate u_rho and v_rho to get kinetic energy budget
#        # (take diag terms, divide by dz, and multiply by u) 
#        u_rho[:,:,j,i] = (unc[:,:,j,i]+unc[:,:,j,i+1])/2
#        v_rho[:,:,j,i] = (vnc[:,:,j,i]+vnc[:,:,j+1,i])/2


'''
dissmag = np.sqrt((udiss_rho**2+vdiss_rho**2))
# divide by dz to get momentum term in m/s^2  (here is hyperdiffusion
# from advective scheme (UP3))
dissmag_dz = dissmag/Hz[:,1:-1,1:-1]
    
vmixmag = np.sqrt((uvmix_rho**2+vvmix_rho**2))
# divide by dz to get momentum term in m/s^2 (here is dissipation)
# d/dz (Akv du/dz)
vmixmag_dz = vmixmag/Hz[:,1:-1,1:-1]

#pgrmag = np.sqrt((upgr_rho**2+vpgr_rho**2))
#pgrmag_dz = pgrmag/Hz[:,1:-1,1:-1]

#advmag = np.sqrt((uadv_rho**2+vadv_rho**2))
#advmag_dz = advmag/Hz[:,1:-1,1:-1]

# calculate velocity magnitude
# get rid of first times step (initial condition)
umag = np.sqrt(u_rho[1:,:,:,:]**2+v_rho[1:,:,:,:]**2)
'''

##################################
# calculate kinetic energy budget
# momentum terms times u
# (diagnostic terms * 1/dz * u)
##################################
#dissbudg = dissmag_dz[-1]*umag
#vmixbudg = vmixmag_dz[-1]*umag
#pgrbudg = pgrmag_dz[-1]*umag
#advbudg = advmag_dz[-1]*umag

diss_dz = udiss/Hz[:,:,:-1]
vmix_dz = uvmix/Hz[:,:,:-1]
#pgr_dz = upgr/Hz[:,:,:-1]
#adv_dz = uadv/Hz[:,:,:-1]

'''
dissbudg_u = udiss*unc[1:,:,:,:]
vmixbudg_u = uvmix*unc[1:,:,:,:]
pgrbudg_u = upgr*unc[1:,:,:,:]
advbudg_u = uadv*unc[1:,:,:,:]
'''
dissbudg_u = udiss*unc
vmixbudg_u = uvmix*unc
#pgrbudg_u = upgr*unc
#advbudg_u = uadv*unc

# 10 m
#dissbudg_x = np.nanmean(dissbudg_u[-1,:,20:110,:],axis=(0,1))
#vmixbudg_x = np.nanmean(vmixbudg_u[-1,:,20:110,:],axis=(0,1)) 
#pgrbudg_x = np.nanmean(pgrbudg_u[-1,:,20:110,:],axis=(0,1))
#advbudg_x = np.nanmean(advbudg_u[-1,:,20:110,:],axis=(0,1))

# 3 m
#dissbudg_x = np.nanmean(dissbudg_u[-1,:,106:406,:],axis=(0,1))
#vmixbudg_x = np.nanmean(vmixbudg_u[-1,:,106:406,:],axis=(0,1))
#pgrbudg_x = np.nanmean(pgrbudg_u[-1,:,106:406,:],axis=(0,1))
#advbudg_x = np.nanmean(advbudg_u[-1,:,106:406,:],axis=(0,1))

# get energy at last time step as a function of x
#dissbudg_x = np.nanmean(dissbudg[-1],axis=(0,1))
#vmixbudg_x = np.nanmean(vmixbudg[-1],axis=(0,1))
#pgrbudg_x = np.nanmean(pgrbudg[-1],axis=(0,1))
#advbudg_x = np.nanmean(advbudg[-1],axis=(0,1))

#plt.plot(xslice_m[1:],dissbudg_x)
#plt.plot(xslice_m[1:],vmixbudg_x)
#plt.plot(xslice_m[1:],pgrbudg_x)
#plt.plot(xslice_m[1:],advbudg_x)

##########
# plotting
##########
figw = 14
figh = 15
axsize = 16

# F = 100, dx = 10 m
#vmin0 = 0
#vmax0 = 1E-3
#
#vmin1 = 0.15
#vmax1 = 0.2

# F = 0, dx = 10 m
vmin0 = 0
vmax0 = 4E-3

vmin1 = 0
vmax1 = 0.015

##################
# magnitude
#################
'''
fig1,ax1 = plt.subplots(2,1,figsize=[figw,figh],constrained_layout=True)
dissplot1 = ax1.flat[0].pcolor(xslice_m[1:-1],zr[:,cent,1:-1],dissmag_dz[-1,:,cent,:],cmap='rainbow',vmin=vmin0,vmax=np.nanmax(dissmag_dz[-1,:,cent,:]))
vmixplot1 = ax1.flat[1].pcolor(xslice_m[1:-1],zr[:,cent,1:-1],vmixmag_dz[-1,:,cent,:],cmap='rainbow',vmin=vmin1,vmax=np.nanmax(vmixmag_dz[-1,:,cent,:]))
#dissplot = ax.flat[0].pcolor(xslice_m[1:-1],zr[:,cent,1:-1],dissmag_dz[-1,:,cent,:],cmap='rainbow')
#vmixplot = ax.flat[1].pcolor(xslice_m[1:-1],zr[:,cent,1:-1],vmixmag_dz[-1,:,cent,:],cmap='rainbow')
#dissplot = ax.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],diss_dz[-1,:,cent,:],cmap='bwr',vmin=-0.0006,vmax=0.0006)
#vmixplot = ax.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmix_dz[-1,:,cent,:],cmap='bwr',vmin=-0.035,vmax=0.035)
ax1.flat[0].set_title('hyperdiffusive term in UP3 '+r'[$\nabla_h\kappa\nabla_h^3\mathbf{u_h}$] (m/s$^2$)',fontsize=axsize)
ax1.flat[1].set_title('vertical mixing magnitude '+r'[$\frac{\partial}{\partial z}(Akv\frac{\partial \mathbf{u}}{\partial z})$] (m/s$^2$)',fontsize=axsize)
#ax.flat[1].set_title('vertical mixing '+r'[$\frac{\partial}{\partial z}(Akv\frac{\partial \mathbf{u}}{\partial z})$] (m/s$^2$)',fontsize=axsize)

ax1.flat[1].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax1.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax1.flat[1].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax1 = fig1.colorbar(dissplot1,ax=ax1.flat[0])
cb1ax1 = fig1.colorbar(vmixplot1,ax=ax1.flat[1])
cb0ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax1.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax1.flat[1].tick_params(axis='both',which='major',labelsize=axsize)
'''

##################
# raw values
#################
fig2,ax2 = plt.subplots(2,1,figsize=[figw,figh],constrained_layout=True)
# 10 m
#dissplot2 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],diss_dz[-1,:,cent,:],cmap='bwr',vmin=-0.0006,vmax=0.0006)
#vmixplot2 = ax2.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmix_dz[-1,:,cent,:],cmap='bwr',vmin=-0.035,vmax=0.035)
#dissplot2 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],diss_dz[-1,:,cent,:],cmap='bwr')
#vmixplot2 = ax2.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmix_dz[-1,:,cent,:],cmap='bwr')
#dissplot2 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],diss_dz[:,cent,:],cmap='bwr',vmin=-2e-4,vmax=2e-4)
#vmixplot2 = ax2.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmix_dz[:,cent,:],cmap='bwr',vmin=-2e-4,vmax=2e-4)
#dissplot2 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],diss_dz[:,cent,:],cmap='bwr',vmin=-0.0006,vmax=0.0006)
#vmixplot2 = ax2.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmix_dz[:,cent,:],cmap='bwr',vmin=-0.035,vmax=0.035)
dissplot2 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],diss_dz[:,cent,:],cmap='bwr',vmin=-0.0006,vmax=0.0006)
vmixplot2 = ax2.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmix_dz[:,cent,:],cmap='bwr',vmin=-0.0001,vmax=0.0001)
ax2.flat[0].set_title('hyperdiffusive term in UP3 '+r'[$\nabla_h\kappa\nabla_h^3\mathbf{u_h}$] (m/s$^2$)',fontsize=axsize)
ax2.flat[1].set_title('vertical mixing '+r'[$\frac{\partial}{\partial z}(A_v\frac{\partial \mathbf{u}}{\partial z})$] (m/s$^2$)',fontsize=axsize)

ax2.flat[1].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax2.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax2.flat[1].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax2 = fig2.colorbar(dissplot2,ax=ax2.flat[0])
cb1ax2 = fig2.colorbar(vmixplot2,ax=ax2.flat[1])
cb0ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax2.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax2.flat[1].tick_params(axis='both',which='major',labelsize=axsize)

#################
# kinetic energy
#################
fig3,ax3 = plt.subplots(2,1,figsize=[figw,figh],constrained_layout=True)
# 10 m
#dissplot3 = ax3.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],dissbudg_u[-1,:,cent,:],cmap='bwr',vmin=-5E-5,vmax=5E-5)
#dissplot3 = ax3.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],dissbudg_u[-1,:,cent,:],cmap='bwr')
#vmixplot3 = ax3.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmixbudg_u[-1,:,cent,:],cmap='rainbow',vmin=0,vmax=np.nanmax(vmixbudg_u[-1,:,cent,:]))
dissplot3 = ax3.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],dissbudg_u[:,cent,:],cmap='bwr',vmin=-2e-4,vmax=2e-4)
vmixplot3 = ax3.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmixbudg_u[:,cent,:],cmap='bwr',vmin=-2e-4,vmax=2e-4)
ax3.flat[0].set_title('hyperdiffusive term in UP3 '+r'[$\kappa(\nabla_h^2\mathbf{u_h})^2]$ (m$^2$/s$^3$)',fontsize=axsize)
ax3.flat[1].set_title('vertical mixing '+r'[Akv$(\frac{\partial \mathbf{u}}{\partial z})^2$] (m$^2$/s$^3$)',fontsize=axsize)
#ax.flat[1].set_title('vertical mixing '+r'[$\frac{\partial}{\partial z}(Akv\frac{\partial \mathbf{u}}{\partial z})$] (m/s$^2$)',fontsize=axsize)

ax3.flat[1].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax3.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax3.flat[1].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax3 = fig3.colorbar(dissplot3,ax=ax3.flat[0])
cb1ax3 = fig3.colorbar(vmixplot3,ax=ax3.flat[1])
cb0ax3.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax3.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax3.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax3.flat[1].tick_params(axis='both',which='major',labelsize=axsize)


#################
# Akv
#################
fig4,ax4 = plt.subplots(1,1,figsize=[figw,figh],constrained_layout=True)
# 10 m
#dissplot3 = ax4.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],dissbudg_u[-1,:,cent,:],cmap='bwr',vmin=-5E-5,vmax=5E-5)
#dissplot3 = ax4.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],dissbudg_u[-1,:,cent,:],cmap='bwr')
#vmixplot3 = ax4.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],vmixbudg_u[-1,:,cent,:],cmap='rainbow',vmin=0,vmax=np.nanmax(vmixbudg_u[-1,:,cent,:]))
akvplt = ax4.pcolor(xslice_m,zw[:,cent,:],akv[:,cent,:],cmap='rainbow',vmin=0,vmax=0.025)
ax4.set_title(r'$A_v$',fontsize=axsize)
#ax.[1].set_title('vertical mixing '+r'[$\frac{\partial}{\partial z}(Akv\frac{\partial \mathbf{u}}{\partial z})$] (m/s$^2$)',fontsize=axsize)

ax4.set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax4.set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.[1],format='%.0e')
cb0ax4 = fig4.colorbar(akvplt,ax=ax4)
cb0ax4.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax4.tick_params(axis='both',which='major',labelsize=axsize)


