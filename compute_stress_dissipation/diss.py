###################
# compute UP3 vs Akv (from gray mix)
# vs vertical advection in implicit step (vmix minus Wimp)
# vs Wimp from plume
###################
import os
import sys
sys.path.append('/data/project3/minnaho/global')
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import ROMS_depths as depths
#plt.ion()

###########
# test cases
###########

#datanc = Dataset('pipe_3m_dt1_hydro_F0_debug_1ts_closed_diag_rnd_his_his.20000101000000.nc','r')
#diagnc = Dataset('pipe_3m_dt1_hydro_F0_debug_1ts_closed_diag_rnd_his_dia.20000101015500.nc','r')
#gridnc = Dataset('pipe_3m_dt1_hydro_F0_debug_1ts_closed_diag_rnd_his_grd.nc','r')
#randnc = Dataset('pipe_3m_dt1_hydro_F0_debug_1ts_closed_diag_rnd_his_rnd.20000101015500.nc','r')

datanc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_his.20000101060000.nc','r')
diagnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_dia.20000101120000.nc','r')
gridnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_grd.nc','r')
randnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_rnd.20000101061500.nc','r')

#datanc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_noplume_his.20000101120000.nc','r')
#diagnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_noplume_dia.20000101120000.nc','r')
#gridnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_noplume_grd.nc','r')
#randnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_noplume_rnd.20000101061500.nc','r')

#savename = 'F1_noplume'
savename = 'typicalsummer_F1'

dx = 3.
#cent = int(128/2)
cent = int(512/2)


#################################################################
# begin reading and calculation
#################################################################

xslice_m = np.array(gridnc.variables['x_rho'])[cent]

udiss = np.array(diagnc.variables['u_dis'])[-1,:,:,:]
vdiss = np.array(diagnc.variables['v_dis'])[-1,:,:,:]

uvmix = np.array(diagnc.variables['u_vmx'])[-1,:,:,:]
vvmix = np.array(diagnc.variables['v_vmx'])[-1,:,:,:]

upgr = np.array(diagnc.variables['u_pgr'])[-1,:,:,:]
vpgr = np.array(diagnc.variables['v_pgr'])[-1,:,:,:]

uadv = np.array(diagnc.variables['u_adv'])[-1,:,:,:]
vadv = np.array(diagnc.variables['v_adv'])[-1,:,:,:]

unc = np.array(datanc.variables['u'])[-1,:,:,:]
vnc = np.array(datanc.variables['v'])[-1,:,:,:]

akv = np.array(datanc.variables['Akv'])[-1,:,:,:]
akvrho = 0.5*(akv[1:]+akv[:-1])
akvw = 0.5*(akvrho[1:]+akvrho[:-1])

akvu = 0.5*(akvw[:,:,1:]+akvw[:,:,:-1])

Wimpnc = np.array(randnc.variables['Wimp'])[-1,:,:,:]

pmnc = np.array(gridnc.variables['pm'])
pnnc = np.array(gridnc.variables['pn'])

wimp = Wimpnc*pmnc*pnnc # turn m^3/s into m/s

udiss[udiss>1E10] = np.nan
vdiss[vdiss>1E10] = np.nan

uvmix[uvmix>1E10] = np.nan
vvmix[vvmix>1E10] = np.nan

upgr[upgr>1E10] = np.nan
vpgr[vpgr>1E10] = np.nan

uadv[uadv>1E10] = np.nan
vadv[vadv>1E10] = np.nan

unc[unc>1E10] = np.nan
vnc[vnc>1E10] = np.nan


# depths and Hz at last time step
zw = depths.get_zw_zeta_tind(datanc,gridnc,-1) # depths at w points
zr = depths.get_zr_zeta_tind(datanc,gridnc,-1) # depths at rho points
Hz = np.diff(zw,axis=0) # height of dz
Hzu = 0.5*(Hz[:,:,1:]+Hz[:,:,:-1]) # get Hz at u points
Hzu_dz = 0.5*(Hzu[1:]+Hzu[:-1])  # get Hzu at w points to match du

# dudz is at w points now [1:-1] (-2 indices)
dudz = (unc[1:]-unc[:-1])/Hzu_dz

# put wimp at u points
wimpu = 0.5*(wimp[:,:,1:]+wimp[:,:,:-1])
# put wimp at rho points
wimpurho = 0.5*(wimpu[1:]+wimpu[:-1])
# put back at w points (-2 indices) to match du/dz
wimpuw = 0.5*(wimpurho[1:]+wimpurho[:-1])

wimpdudz = wimpuw*dudz

# diss
# divide by dz to get momentum term in m/s^2  (here is hyperdiffusion
# from horizontal advective scheme (UP3))
udiss_dz = udiss/Hzu

# adv
# divide by dz to get momentum term in m/s^2  (here is 4th order
# centered advection term
# from horizontal advective scheme (UP3))
uadv_dz = uadv/Hzu

# UP3
uup3 = uadv_dz+udiss_dz
    
# vmix
# divide by dz to get momentum term in m/s^2 (here is dissipation)
# d/dz (Akv du/dz)
uvmix_dz = uvmix/Hzu

# vmix at w points (-2 indices to match wimpdudz)
uvmixw_dz = 0.5*(uvmix_dz[1:]+uvmix_dz[:-1])

# pgr
# divide by dz to get momentum term in m/s^2 (here is pressure
# gradient term)
# -del_h(phi)
upgr_dz = upgr/Hzu


##########
# plotting
##########
figw = 14
figh = 15
axsize = 16

####################
# plot diffusion, advection, and 
# UP3 = diss + adv
####################
fig1,ax1 = plt.subplots(3,1,figsize=[figw,figh],constrained_layout=True,sharex=True)
dissplot1 = ax1.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],udiss_dz[:,cent,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
advplot1 = ax1.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],uadv_dz[:,cent,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
up3plot1 = ax1.flat[2].pcolor(xslice_m[:-1],zr[:,cent,:-1],uup3[:,cent,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
ax1.flat[0].set_title('hyperdiffusive term in UP3 '+r'[$\nabla_h\kappa\nabla_h^3\mathbf{u_h}$] (m/s$^2$)',fontsize=axsize)
ax1.flat[1].set_title('centered advection term in UP3 '+r'[$\nabla_h(\mathbf{u_h}^2)$] (m/s$^2$)',fontsize=axsize)
ax1.flat[2].set_title('UP3 = advection + hyperdiffusion',fontsize=axsize)

ax1.flat[2].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax1.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax1.flat[1].set_ylabel('depth (m)',fontsize=axsize)
ax1.flat[2].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax1 = fig1.colorbar(dissplot1,ax=ax1.flat[0],format='%.0e')
cb1ax1 = fig1.colorbar(advplot1,ax=ax1.flat[1],format='%.0e')
cb2ax1 = fig1.colorbar(up3plot1,ax=ax1.flat[2],format='%.0e')
cb0ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb2ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax1.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax1.flat[1].tick_params(axis='both',which='major',labelsize=axsize)
ax1.flat[2].tick_params(axis='both',which='major',labelsize=axsize)

fig1.savefig(savename+'_up3.png')


####################
# plot vmx (d/dz(Akv(du/dz))),
# Wimp, vmx-Wimp (vertical advection in implicit step)
####################
fig2,ax2 = plt.subplots(3,1,figsize=[figw,figh],constrained_layout=True,sharex=True)
vmixplot1 = ax2.flat[0].pcolor(xslice_m[:-1],zw[1:-1,cent,:-1],uvmixw_dz[:,cent,:],cmap='bwr',vmin=-5E-5,vmax=5E-5)
wimpplot1 = ax2.flat[1].pcolor(xslice_m[:-1],zw[1:-1,cent,:-1],wimpdudz[:,cent,:],cmap='bwr',vmin=-1E-19,vmax=1E-19)
diffplot1 = ax2.flat[2].pcolor(xslice_m[:-1],zw[1:-1,cent,:-1],uvmixw_dz[:,cent,:]-wimpdudz[:,cent,:],cmap='bwr',vmin=-5E-5,vmax=5E-5)
ax2.flat[0].set_title('vertical mixing '+r'$\frac{\partial}{\partial z}$(Akv$\frac{\partial u}{\partial z})$ (m/s$^2$)',fontsize=axsize)
ax2.flat[1].set_title('implicit omega '+r'$w\frac{\partial u}{\partial z}$ (m/s$^2$)',fontsize=axsize)
ax2.flat[2].set_title('vertical mixing - implicit omega',fontsize=axsize)

ax2.flat[2].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax2.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax2.flat[1].set_ylabel('depth (m)',fontsize=axsize)
ax2.flat[2].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax2 = fig2.colorbar(vmixplot1,ax=ax2.flat[0],format='%.0e')
cb1ax2 = fig2.colorbar(wimpplot1,ax=ax2.flat[1],format='%.0e')
cb2ax2 = fig2.colorbar(diffplot1,ax=ax2.flat[2],format='%.0e')
cb0ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb2ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax2.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax2.flat[1].tick_params(axis='both',which='major',labelsize=axsize)
ax2.flat[2].tick_params(axis='both',which='major',labelsize=axsize)

fig2.savefig(savename+'_vmix_Wimp.png')

fig3,ax3 = plt.subplots(1,1,figsize=[figw,figh*(1./3)],constrained_layout=True)
akvplot = ax3.pcolor(xslice_m[:-1],zw[1:-1,cent,:-1],akvu[:,cent,:],cmap='rainbow',vmin=0,vmax=1E-3)
ax3.set_title('Akv '+r'(m$^2$/s)',fontsize=axsize)
ax3.set_ylabel('depth (m)',fontsize=axsize)
ax3.set_xlabel('cross pipe distance (m)',fontsize=axsize)
cb0ax3 = fig3.colorbar(akvplot,ax=ax3)
ax3.tick_params(axis='both',which='major',labelsize=axsize)
cb0ax3.ax.tick_params(axis='both',which='major',labelsize=axsize)


fig3.savefig(savename+'_akv_uw.png')


# calculate dissipation
fig4,ax4 = plt.subplots(3,1,figsize=[figw,figh],constrained_layout=True,sharex=True)
dissplot1 = ax4.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],udiss_dz[:,cent,:]*unc[:,cent,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
advplot1 = ax4.flat[1].pcolor(xslice_m[:-1],zr[:,cent,:-1],uadv_dz[:,cent,:]*unc[:,cent,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
up3plot1 = ax4.flat[2].pcolor(xslice_m[:-1],zr[:,cent,:-1],uup3[:,cent,:]*unc[:,cent,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
ax4.flat[0].set_title('hyperdiffusive term in UP3 * u',fontsize=axsize)
ax4.flat[1].set_title('centered advection term in UP3 * u',fontsize=axsize)
ax4.flat[2].set_title('dissipation from UP3 '+r'$\varepsilon$ (m$^2$/s$^3$)',fontsize=axsize)

ax4.flat[2].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax4.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax4.flat[1].set_ylabel('depth (m)',fontsize=axsize)
ax4.flat[2].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax4 = fig4.colorbar(dissplot1,ax=ax4.flat[0],format='%.0e')
cb1ax4 = fig4.colorbar(advplot1,ax=ax4.flat[1],format='%.0e')
cb2ax4 = fig4.colorbar(up3plot1,ax=ax4.flat[2],format='%.0e')
cb0ax4.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax4.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb2ax4.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax4.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax4.flat[1].tick_params(axis='both',which='major',labelsize=axsize)
ax4.flat[2].tick_params(axis='both',which='major',labelsize=axsize)

fig4.savefig(savename+'_disspation.png')

