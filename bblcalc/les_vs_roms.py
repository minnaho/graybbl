import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import ROMS_depths as depths

plt.ion()

###########
# LES data
###########
# PGF = u_g*f = 6x10^-6
datastr = dict(np.load('flat_BBL_f_strat.npz',allow_pickle=True))
datawif = dict(np.load('flat_BBL_f.npz',allow_pickle=True))

# time start to end to average over
# 180 = last 10 hours
# 156 = last 12 hours
timeslice = 156

str_dz    = datastr['dz']
str_dx    = datastr['dx']
str_zplt  = str_dz*np.arange(1,datastr['nz']+str_dz)
str_xplt  = str_dx*np.arange(1,datastr['nx']+str_dx)
str_tplt  = np.arange(timeslice+1,datastr['u_prof'].shape[0]+1)*300/3600 # time in minutes
str_tplt_full  = np.arange(1,datastr['u_prof'].shape[0]+1)*300/3600 # time in minutes
str_avg_u = datastr['avg_u']
str_uprof = np.nanmean(datastr['u_prof'][timeslice:],axis=0)
str_upwp = np.nanmean(datastr['upwp_prof'][timeslice:],axis=0)

# unresolved stresses
#str_unresolved = np.sqrt(np.nanmean(datastr['txz_prof'][timeslice:],axis=0)**2+np.nanmean(datastr['tyz_prof'][timeslice:],axis=0)**2)
str_unresolved = np.nanmean(datastr['txz_prof'][timeslice:],axis=0)

###########
# ROMS data
###########

#datanc = Dataset('obc10m/fb_stratified_gray_obc_10m_his_concat.nc','r')
#datagrd = Dataset('obc10m/fb_stratified_gray_obc_10m_grd.nc','r')
#datarnd = Dataset('obc10m/fb_stratified_gray_obc_10m_rnd_concat.nc','r')


#datanc =  Dataset('ewp10m/fb_strat_gray_ewp_nso_10m_his_concat.nc','r') 
#datagrd = Dataset('ewp10m/fb_strat_gray_ewp_nso_10m_grd.nc','r')
#datarnd = Dataset('ewp10m/fb_strat_gray_ewp_nso_10m_rnd_concat.nc','r')

#datanc =  Dataset('lesinit/fb_nh_strat_gray_ewp_nso_lesinit_10m_his_concat.nc','r') 
#datagrd = Dataset('lesinit/fb_nh_strat_gray_ewp_nso_lesinit_10m_grd.nc','r')
#datarnd = Dataset('lesinit/fb_nh_strat_gray_ewp_nso_lesinit_10m_rnd_concat.nc','r')

datanc =  Dataset('lesinit/fb_nh_strat_ewp_nso_lesinit_10m_his_concat.nc','r') 
datagrd = Dataset('lesinit/fb_nh_strat_ewp_nso_lesinit_10m_grd.nc','r')
datarnd = Dataset('lesinit/fb_nh_strat_ewp_nso_lesinit_10m_rnd_concat.nc','r')

roms_x = datagrd.variables['x_rho']

# u'w' calculation
roms_unc_full = np.squeeze(datanc.variables['u'])
roms_unc = np.squeeze(datanc.variables['u'][timeslice+1:])
roms_wnc = np.squeeze(datanc.variables['w'][timeslice+1:])
roms_vnc = np.squeeze(datanc.variables['v'][timeslice+1:])
roms_zetanc = np.squeeze(datanc.variables['zeta'][timeslice+1:])

roms_tplt = np.arange(timeslice+1,roms_unc.shape[0]+timeslice+1)*300/3600
roms_tplt_full = np.arange(0,roms_unc.shape[0]+timeslice+1)*300/3600
roms_uhovplt = np.nanmean(roms_unc,axis=(2,3))
roms_uhovplt_full = np.nanmean(roms_unc_full,axis=(2,3))

roms_wavgtime = np.nanmean(roms_wnc,axis=0)
roms_uavgtime = np.nanmean(roms_unc,axis=0)
roms_zetatime = np.nanmean(roms_zetanc,axis=0)

zr = depths.get_zr_zeta(datanc,datagrd,roms_zetatime)
zr_xymean = np.nanmean(zr[:,10:-10,10:-10],axis=(1,2))

roms_uavg_xymean = np.nanmean(roms_uavgtime[:,10:-10,10:-10],axis=(1,2))

roms_up = roms_unc-roms_uavgtime
roms_wp = roms_wnc-roms_wavgtime

roms_wp_u = 0.5*(roms_wp[:,:,:,1:]+roms_wp[:,:,:,:-1])

roms_upwp = roms_up*roms_wp_u
roms_upwp_avg = np.nanmean(roms_upwp[:,:,10:-10,10:-10],axis=(0,2,3))

# bottom drag
rd = np.squeeze(datarnd.variables['rd_out'][timeslice:])
rd_uv1 = 0.5*(rd[:,:,1:]+rd[:,:,:-1])
rd_uv = 0.5*(rd_uv1[:,1:,:]+rd_uv1[:,:-1,:])

#tau_roms = np.sqrt((rd_uv*0.5*(roms_unc[:,0,1:,:]+roms_unc[:,0,:-1,:]))**2+(rd_uv*0.5*(roms_vnc[:,0,:,1:]+roms_vnc[:,0,:,:-1]))**2)
tau_roms = rd_uv1*0.5*(roms_unc[:,0,:,:]+roms_unc[:,0,:,:])
tau_mean = np.ones((roms_unc.shape[1]-1))*np.nanmean(tau_roms,axis=(0,1,2))
tau_mean[1:] = 0

roms_rho = np.squeeze(datanc.variables['rho'][timeslice+1:])

roms_rhoavgtime = np.nanmean(roms_rho,axis=0)
roms_vavgtime = np.nanmean(roms_vnc,axis=0)

roms_rho_uv = 0.5*(roms_rhoavgtime[:,1:,1:]+roms_rhoavgtime[:,:-1,:-1])+1000

vonKar = 0.41
Zob = 1.0E-2

Cd = (vonKar**2)/(np.log((zr_xymean[0]+48)/Zob))**2

roms_uavgtime_v = 0.5*(roms_uavgtime[:,1:,:]+roms_uavgtime[:,:-1,:])
roms_vavgtime_u = 0.5*(roms_vavgtime[:,:,1:]+roms_vavgtime[:,:,:-1])

roms_ub = np.sqrt(roms_uavgtime_v**2+roms_vavgtime_u)

tau = Cd*roms_ub**2

tau_stress = 0.5*(tau[1:]+tau[:-1])

#tau_xymean = np.nanmean(tau_stress,axis=(1,2))

# unresolved stress calculation (kv*du/dz)
akvnc = np.squeeze(datanc.variables['Akv'][timeslice:])
akvnc_avg_xymean = np.nanmean(akvnc[:,:,10:-10,10:-10],axis=(0,2,3))

akvnc_avg_u1 = 0.5*(akvnc_avg_xymean[1:]+akvnc_avg_xymean[:-1])
akvnc_avg_u = 0.5*(akvnc_avg_u1[1:]+akvnc_avg_u1[:-1])

zr_u = 0.5*(zr[:,:,1:]+zr[:,:,:-1])

# multiply dz by -1 because shouldn't have negative distances
dudz = np.nanmean(roms_unc[:,:-1,:,:]-roms_unc[:,1:,:,:],axis=0)/((zr_u[:-1,:,:]-zr_u[1:,:,:])*-1)

dudz_avg = np.nanmean(dudz,axis=(1,2))

kvdudz = akvnc_avg_u*dudz_avg

zr_dudz = 0.5*(zr_xymean[1:]+zr_xymean[:-1])

###########
# plotting
###########

#fig1,ax1 = plt.subplots(1,3,figsize=[15,6])
#ax1.flat[0].plot(wif_uprof,wif_zplt)
#ax1.flat[0].set_xlabel('u')
#ax1.flat[1].plot(wif_dpdx,wif_zplt)
#ax1.flat[1].set_xlabel('dpdx')
#ax1.flat[2].plot(wif_avgdpdx,wif_zplt)
#ax1.flat[2].set_xlabel('avg_dpdx')

# average u hovmoller
#fig2,ax2 = plt.subplots(2,1,figsize=[15,10])
## last time step of u (averged over last 5 minutes)
#avguplt2 = ax2.flat[0].pcolormesh(str_xplt,str_zplt,str_avg_u[-1],cmap='rainbow',vmin=0,vmax=0.061)
#cb2ax0 = fig2.colorbar(avguplt2,ax=ax2.flat[0])
##uprofplt2 = ax2.flat[1].pcolormesh(str_tplt,str_zplt,str_uprof.transpose(),cmap='rainbow',vmin=0,vmax=0.061)
##cb2ax1 = fig2.colorbar(uprofplt2,ax=ax2.flat[1])
#avgupltzoom2 = ax2.flat[1].pcolormesh(str_xplt[64:128],str_zplt[:31],str_avg_u[-1,64:128,:31].transpose(),cmap='rainbow',vmin=0,vmax=0.061)
#cb2ax1 = fig2.colorbar(avguplt2,ax=ax2.flat[1])
#
#fig3,ax3 = plt.subplots(2,1,figsize=[15,10])
#avguplt3 = ax3.flat[0].pcolormesh(roms_x[65],zr[:,65,:]+48,roms_uavgtime[:,65,:],cmap='rainbow',vmin=0,vmax=0.061)
#cb3ax0 = fig3.colorbar(avguplt3,ax=ax3.flat[0])
#avgupltzoom3 = ax3.flat[1].pcolormesh(roms_x[65,86:172],zr[:12,65,86:172]+48,roms_uavgtime[:12,65,86:172],cmap='rainbow',vmin=0,vmax=0.061)
#cb3ax1 = fig3.colorbar(avgupltzoom3,ax=ax3.flat[1])

# u hovmoller
fig1,ax1 = plt.subplots(2,1,sharey=True,sharex=True,figsize=[15,10])
les_uhov = ax1.flat[0].pcolormesh(str_tplt,str_zplt,np.transpose(datastr['u_prof'][timeslice:]),cmap='rainbow',vmin=0.03,vmax=0.061)
cb1ax0 = fig1.colorbar(les_uhov,ax=ax1.flat[0])
roms_uhov = ax1.flat[1].pcolormesh(roms_tplt,zr_xymean+48,np.transpose(roms_uhovplt),cmap='rainbow',vmin=0.03,vmax=0.061)
cb1ax1 = fig1.colorbar(roms_uhov,ax=ax1.flat[1])
ax1.flat[1].set_xlabel('hours')
ax1.flat[0].set_ylabel('depth (m)')
ax1.flat[1].set_ylabel('depth (m)')
ax1.flat[0].set_title('LES')
ax1.flat[1].set_title('ROMS')

fig2,ax2 = plt.subplots(2,1,sharey=True,sharex=True,figsize=[15,10])
les_uhov = ax2.flat[0].pcolormesh(str_tplt_full,str_zplt,np.transpose(datastr['u_prof']),cmap='rainbow',vmin=0.03,vmax=0.061)
cb2ax0 = fig2.colorbar(les_uhov,ax=ax2.flat[0])
roms_uhov = ax2.flat[1].pcolormesh(roms_tplt_full,zr_xymean+48,np.transpose(roms_uhovplt_full),cmap='rainbow',vmin=0.03,vmax=0.061)
cb2ax1 = fig2.colorbar(roms_uhov,ax=ax2.flat[1])
ax2.flat[1].set_xlabel('hours')
ax2.flat[0].set_ylabel('depth (m)')
ax2.flat[1].set_ylabel('depth (m)')
ax2.flat[0].set_title('LES')
ax2.flat[1].set_title('ROMS')


# u profile compare
fig4,ax4 = plt.subplots(1,1,sharey=True,sharex=True,figsize=[6,6])
ax4.plot(str_uprof,str_zplt,label='LES')
ax4.set_xlabel('u (m/s)')
ax4.set_ylabel('depth (m)')

ax4.plot(roms_uavg_xymean,zr_xymean+48,label='ROMS')
ax4.legend()

# u'w'
# "resolved stresses"
fig5,ax5 = plt.subplots(1,1,sharex=True,sharey=True,figsize=[6,6])
ax5.plot(str_upwp-str_unresolved,str_zplt,label='LES')
ax5.set_xlabel(r'$\overline{u^{\prime}w^{\prime}}$')
ax5.set_ylabel('depth (m)')

ax5.plot(roms_upwp_avg,zr_xymean+48,label='ROMS')
ax5.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
ax5.legend()

# unresolved stresses
fig6,ax6 = plt.subplots(1,1,sharey=True,sharex=True,figsize=[6,6])
ax6.plot(str_unresolved,str_zplt,label='LES '+r'$\tau_{xz}$')
ax6.set_xlabel('sub-grid scale stresses')
ax6.set_ylabel('depth (m)')

ax6.plot(kvdudz-tau_mean,zr_dudz+48,label='ROMS '+r'$kv\frac{\partial u}{\partial z}-\tau_x$')
ax6.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
ax6.legend()

# resolved + unresolved stresses
fig7,ax7 = plt.subplots(1,1,sharey=True,sharex=True,figsize=[6,6])
ax7.plot(str_upwp,str_zplt,label='LES')
ax7.set_xlabel('resolved + unresolved stresses')
ax7.set_ylabel('depth (m)')

#ax7.flat[1].set_xlabel('resolved + unresolved stresses')
#ax7.flat[1].set_title('ROMS')
ax7.plot(0.5*(roms_upwp_avg[1:]+roms_upwp_avg[:-1])+kvdudz-tau_mean,zr_dudz+48,label='ROMS')
ax7.legend()
ax7.ticklabel_format(axis='x',style='sci',scilimits=(0,0))

