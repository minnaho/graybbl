###################
# compute bar(w'^2), bar(T'^2) (variance,numpy.var)
# dissipation (using diag?), bar(u) (time)
# covariance to get bar(u'w') (np.cov(x,y))
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

#datanc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_his.20000101060000.nc','r')
#diagnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_dia.20000101120000.nc','r')
#gridnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F1_typicalsummer_grd.nc','r')

datanc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F01_typicalsummer_his.20000101060000.nc','r')
diagnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F01_typicalsummer_dia.20000101120000.nc','r')
gridnc = Dataset('/data/project9/minnaho/graybbl/movies/pipe_3m_dt1_nh_F01_typicalsummer_grd.nc','r')

dx = 3.
cent = int(512/2)


#################################################################
# begin reading and calculation
#################################################################

xslice_m = np.array(gridnc.variables['x_rho'])[cent]

udiss = np.array(diagnc.variables['u_dis'])[-1,:,:,:]

unc = np.array(datanc.variables['u'])[:,:,:,:]
#vnc = np.array(datanc.variables['v'])[:,:,:,:]
wnc = np.array(datanc.variables['w'])[:,:,:,:]

tempnc = np.array(datanc.variables['temp'])[:,:,:,:]

#akv = np.array(datanc.variables['Akv'])[-1,:,:,:]

udiss[udiss>1E10] = np.nan

unc[unc>1E10] = np.nan
#vnc[vnc>1E10] = np.nan
wnc[wnc>1E10] = np.nan

tempnc[tempnc>1E10] = np.nan

'''
vonKar = 0.41
Zob = 1E-2
C_d = (vonKar/np.log(1+0.5*Hz[0,:,:]/Zob))**2
ubot = unc[0,:,:]
rho = np.array(datanc.variables['rho'])[-1,0,:,:] # bottom density
rhou = 0.5*(rho[:,1:]+rho[:,:-1])
C_du = 0.5*(C_d[:,1:]+C_d[:,:-1])

akv_dudz_dz = np.empty((unc.shape[0],unc.shape[1],unc.shape[2]))*np.nan
akv_u = 0.5*(akv[:,:,:-1]+akv[:,:,1:])
Hz_u = 0.5*(Hz[:,:,:-1]+Hz[:,:,1:])
akv_dudz_dz[0,:,:] = rhou*C_du*ubot # bottom stress
for k in range(akv_dudz_dz.shape[0]-2):
    akv_dudz_dz[k+1,:,:] = akv_u[k+1,:,:]*(unc[k+1,:,:]-unc[k,:,:])/Hz_u[k+1,:,:]-akv_dudz_dz[k,:,:]
'''


# depths and Hz at last time step
zw = depths.get_zw_zeta_tind(datanc,gridnc,-1) # depths at w points
zr = depths.get_zr_zeta_tind(datanc,gridnc,-1) # depths at rho points
Hz = np.diff(zw,axis=0) # height of dz
Hzu = 0.5*(Hz[:,:,1:]+Hz[:,:,:-1]) # get Hz at u points

# diss
# divide by dz to get momentum term in m/s^2  (here is hyperdiffusion
# from horizontal advective scheme (UP3))
udiss_dz = udiss/Hzu

# averages
ubar = np.nanmean(unc,axis=0)
tempbar = np.nanmean(tempnc,axis=0)

# variances
wvar = np.var(wnc,axis=0)
tempvar = np.var(tempnc,axis=0)
wu = 0.5*(wnc[:,:,:,1:]+wnc[:,:,:,:-1])

uwcov = np.empty((unc.shape[1],unc.shape[3]))*np.nan

for k_i in range(unc.shape[1]):
    print(str(k_i))
    for i_i in range(unc.shape[3]):
        uwcov[k_i,i_i] = np.cov(unc[:,k_i,cent,i_i],wu[:,k_i,cent,i_i])[0,1]

##########
# plotting
##########
figw = 14
figh = 15
axsize = 16

####################
# plot variance
# bar(w'^2)
####################
fig1,ax1 = plt.subplots(3,1,figsize=[figw,figh],constrained_layout=True,sharex=True)
wvarplot1 = ax1.flat[0].pcolor(xslice_m,zr[:,cent],wvar[:,cent,:],cmap='rainbow',vmin=0,vmax=5E-5)
tempplot1 = ax1.flat[1].pcolor(xslice_m,zr[:,cent],tempvar[:,cent,:],cmap='rainbow',vmin=0,vmax=0.03)
uwcovplot1 = ax1.flat[2].pcolor(xslice_m[:-1],zr[:,cent,:-1],uwcov[:,:],cmap='bwr',vmin=-1E-4,vmax=1E-4)
ax1.flat[0].set_title(r'$\overline{w^{\prime2}}$',fontsize=axsize)
ax1.flat[1].set_title(r'$\overline{T^{\prime2}}$',fontsize=axsize)
ax1.flat[2].set_title(r'$\overline{u^{\prime}w^{\prime}}$',fontsize=axsize)

ax1.flat[2].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax1.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax1.flat[1].set_ylabel('depth (m)',fontsize=axsize)
ax1.flat[2].set_ylabel('depth (m)',fontsize=axsize)

#fig.colorbar(dissplot,ax=ax.flat[0],format='%.0e')
#fig.colorbar(vmixplot,ax=ax.flat[1],format='%.0e')
cb0ax1 = fig1.colorbar(wvarplot1,ax=ax1.flat[0],format='%.0e')
cb1ax1 = fig1.colorbar(tempplot1,ax=ax1.flat[1])
cb2ax1 = fig1.colorbar(uwcovplot1,ax=ax1.flat[2],format='%.0e')
cb0ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb2ax1.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax1.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax1.flat[1].tick_params(axis='both',which='major',labelsize=axsize)
ax1.flat[2].tick_params(axis='both',which='major',labelsize=axsize)

fig1.savefig('var_F01_typicalsummer.png')

##############################
# plot averages (u bar, T bar)
##############################
fig2,ax2 = plt.subplots(3,1,figsize=[figw,figh],constrained_layout=True,sharex=True)
ubarplot1 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],ubar[:,cent,:],cmap='bwr',vmin=-0.15,vmax=0.15)
#ubarplot1 = ax2.flat[0].pcolor(xslice_m[:-1],zr[:,cent,:-1],ubar[:,cent,:],cmap='bwr',vmin=-0.2,vmax=0.2)
tempplot1 = ax2.flat[1].pcolor(xslice_m,zr[:,cent],tempbar[:,cent,:],cmap='rainbow',vmin=11,vmax=18)
udissplot1 = ax2.flat[2].pcolor(xslice_m,zr[:,cent],udiss_dz[:,cent,:],cmap='bwr',vmin=-2E-6,vmax=2E-6)
ax2.flat[0].set_title(r'$\overline{u}$',fontsize=axsize)
ax2.flat[1].set_title(r'$\overline{T}$',fontsize=axsize)
ax2.flat[2].set_title('hyperdiffusive term in UP3 '+r'[$\nabla_h\kappa\nabla_h^3\mathbf{u_h}$] (m/s$^2$)',fontsize=axsize)

ax2.flat[2].set_xlabel('cross pipe distance (m)',fontsize=axsize)

ax2.flat[0].set_ylabel('depth (m)',fontsize=axsize)
ax2.flat[1].set_ylabel('depth (m)',fontsize=axsize)
ax2.flat[2].set_ylabel('depth (m)',fontsize=axsize)

cb0ax2 = fig2.colorbar(ubarplot1,ax=ax2.flat[0])
cb1ax2 = fig2.colorbar(tempplot1,ax=ax2.flat[1])
cb2ax2 = fig2.colorbar(udissplot1,ax=ax2.flat[2],format='%.0e')
cb0ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb1ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)
cb2ax2.ax.tick_params(axis='both',which='major',labelsize=axsize)

ax2.flat[0].tick_params(axis='both',which='major',labelsize=axsize)
ax2.flat[1].tick_params(axis='both',which='major',labelsize=axsize)
ax2.flat[2].tick_params(axis='both',which='major',labelsize=axsize)

fig2.savefig('bar_F01_typicalsummer.png')

