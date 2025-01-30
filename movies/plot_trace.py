# make vertical velocity (w) comparisons
# along pipe
# mean and RMS (root mean square)
import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
from netCDF4 import Dataset
import ROMS_depths as depths
#import seawater as sw
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import glob

# plume, linearly stratified
#output_path0 = 'pipes_3m_nh_dt1_grayn1_F1_his.20000101000000.nc' 
#output_path0 = 'pipes_3m_nh_dt1_grayn1_F1_his.20000101060000.nc' 
#output_path0 = '/data/project9/minnaho/graybbl/F100_tests/pipes_10m_nh_gray_impnoslip_n1_F100_his.20000101000000.nc' 
#output_path0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_noslip_his.20000101000000.nc' 
#output_path0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_3m_dt1_nh_F100_gray_his.20000101000000.nc' 
output_path0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_gray_max_his.20000101000000.nc' 
#output_path0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_kpp_his.20000101000000.nc' 
# plume, stratified
#output_path1 = 'pipe_3m_dt1_nh_F1_typicalsummer_his.20000101000000.nc' 
#output_path1 = 'pipe_3m_dt1_nh_F1_typicalsummer_his.20000101060000.nc' 

# weak shear
#output_path2 = 'pipe_3m_dt1_nh_F0_F01_his.20000101000000.nc'
#output_path2 = 'pipe_3m_dt1_nh_F0_F01_his.20000101060000.nc'
# strong shear
#output_path3 = 'pipe_3m_dt1_nh_F01_F1_his.20000101000000.nc'
#output_path3 = 'pipe_3m_dt1_nh_F01_F1_his.20000101060000.nc'

#grd0 = 'pipes_3m_nh_dt1_grayn1_F1_grd.nc'
#grd0 = '/data/project9/minnaho/graybbl/F100_tests/pipes_10m_nh_gray_impnoslip_n1_F100_grd.nc'
#grd0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_3m_dt1_nh_F100_gray_grd.nc'
grd0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_gray_max_grd.nc'
#grd0 = '/data/project9/minnaho/graybbl/F100_tests/pipe_10m_nh_F100_kpp_grd.nc'
#grd1 = 'pipe_3m_dt1_nh_F1_typicalsummer_grd.nc'
#grd2 = 'pipe_3m_dt1_nh_F0_F01_grd.nc'
#grd3 = 'pipe_3m_dt1_nh_F01_F1_grd.nc'

# plume, F01, stratified
#output_path0 = 'pipe_3m_dt1_nh_F01_typicalsummer_his.20000101000000.nc'
#output_path0 = 'pipe_3m_dt1_nh_F01_typicalsummer_his.20000101060000.nc'
# plume, F1, stratified
#output_path1 = 'pipe_3m_dt1_nh_F1_typicalsummer_his.20000101000000.nc'
output_path1 = 'pipe_3m_dt1_nh_F1_typicalsummer_his.20000101060000.nc'

# no plume, F01 F1 shear, stratified
#output_path2 = 'pipe_3m_dt1_nh_F01_F1_noplume_his.20000101000000.nc'
output_path2 = 'pipe_3m_dt1_nh_F01_F1_noplume_his.20000101060000.nc'
# plume, F01 F1 shear, stratified
#output_path3 = 'pipe_3m_dt1_nh_F0_F01_his.20000101000000.nc'
output_path3 = 'pipe_3m_dt1_nh_F0_F01_his.20000101060000.nc'

#grd0 = 'pipe_3m_dt1_nh_F01_typicalsummer_grd.nc'
grd1 = 'pipe_3m_dt1_nh_F1_typicalsummer_grd.nc'
grd2 = 'pipe_3m_dt1_nh_F01_F1_noplume_grd.nc'
grd3 = 'pipe_3m_dt1_nh_F01_F1_grd.nc'

#cent = 256
#quart = 256-20
#cent = 64
#quart = 128-5
cent = 64
quart = 64-5


xslice_m = np.array(Dataset(grd0,'r').variables['x_rho'])[cent]

#pipe_ls = 106
#pipe_le = 405
pipe_ls = 19
pipe_le = 109

# pick file (0024)
#output0 = list(sorted(glob.glob(output_path0+'*')))[-1]
#output1 = list(sorted(glob.glob(output_path1+'*')))[-4]
#output2 = list(sorted(glob.glob(output_path2+'*')))[-4]
#output3 = list(sorted(glob.glob(output_path3+'*')))[-1]
#output4 = list(sorted(glob.glob(output_path4+'*')))[-1]

# get tracer
tpas_nc0 = np.array(Dataset(output_path0,'r').variables['trace1'])
tpas_nc1 = np.array(Dataset(output_path1,'r').variables['trace1'])
tpas_nc2 = np.array(Dataset(output_path2,'r').variables['trace1'])
tpas_nc3 = np.array(Dataset(output_path3,'r').variables['trace1'])
#tpas_nc4 = Dataset(output_path4,'r').variables['tpas'][:,:,:,:]

#tpas_avg4 = tpas_nc4[tstep,:,pipe_cl,255:]

#tpas_avg0 = np.mean(np.mean(tpas_nc0[:,:,1:513,:],axis=0),axis=1)
#tpas_avg1 = np.mean(np.mean(tpas_nc1[:,:,1:513,:],axis=0),axis=1)
#tpas_avg2 = np.mean(np.mean(tpas_nc2[:,:,1:513,:],axis=0),axis=1)
#tpas_avg3 = np.mean(np.mean(tpas_nc3[:,:,1:513,:],axis=0),axis=1)
#tpas_avg4 = np.mean(np.mean(tpas_nc4[:,:,1:513,:],axis=0),axis=1)


######################
# plotting
#################

#cma = 'bwr'
#cma = 'seismic'
cma = 'gnuplot2_r'
#plt.ion()
axisfont = 16
axistick = 14
vmit = 0
vmat = 10

#savename0 = './figs_F1/strat_F1_trace_'
#savename2 = './figs_shear/shear_weak_strong_trace_'
savename0 = 'F100__'
savename2 = 'F100_'


for t_i in range(Dataset(output_path0,'r').variables['trace1'].shape[0]):
    fig0,(ax0,ax1) = plt.subplots(2,1,figsize=[12,10])
    fig2,(ax2,ax3) = plt.subplots(2,1,figsize=[12,10])
    zr0 = depths.get_zr_zeta_tind(Dataset(output_path0,'r'),Dataset(grd0,'r'),t_i)
    zr1 = depths.get_zr_zeta_tind(Dataset(output_path1,'r'),Dataset(grd1,'r'),t_i)
    zr2 = depths.get_zr_zeta_tind(Dataset(output_path2,'r'),Dataset(grd2,'r'),t_i)
    zr3 = depths.get_zr_zeta_tind(Dataset(output_path3,'r'),Dataset(grd3,'r'),t_i)

    # subtract quart*3 to get m from pipe
    plot0 = ax0.pcolor(xslice_m[quart:]-((quart+5)*10),zr0[:,cent,quart:]*-1,tpas_nc0[t_i,:,cent,quart:],cmap=cma,vmin=vmit,vmax=vmat)
    plot1 = ax1.pcolor(xslice_m[quart:]-((quart+20)*3),zr1[:,cent,quart:]*-1,tpas_nc1[t_i,:,cent,quart:],cmap=cma,vmin=vmit,vmax=vmat)
    plot2 = ax2.pcolor(xslice_m[quart:]-((quart+20)*3),zr2[:,cent,quart:]*-1,tpas_nc2[t_i,:,cent,quart:],cmap=cma,vmin=vmit,vmax=vmat)
    plot3 = ax3.pcolor(xslice_m[quart:]-((quart+20)*3),zr3[:,cent,quart:]*-1,tpas_nc3[t_i,:,cent,quart:],cmap=cma,vmin=vmit,vmax=vmat)

    ax0.invert_yaxis()
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()

    p0 = ax0.get_position().get_points().flatten()
    p1 = ax1.get_position().get_points().flatten()
    p2 = ax2.get_position().get_points().flatten()
    p3 = ax3.get_position().get_points().flatten()

    cbar_ax0 = fig0.add_axes([p1[2]+.015,p1[1],.01,p0[3]-p1[1]])
    cbar_ax2 = fig2.add_axes([p3[2]+.015,p3[1],.01,p2[3]-p3[1]])

    cbar0 = fig0.colorbar(plot0,cax=cbar_ax0,orientation='vertical')
    cbar2 = fig2.colorbar(plot2,cax=cbar_ax2,orientation='vertical')
    cbar0.set_label('tracer (mmol m$^{-3}$)',fontsize=axisfont)
    cbar2.set_label('tracer (mmol m$^{-3}$)',fontsize=axisfont)
    cbar0.ax.tick_params(labelsize=axistick)
    cbar2.ax.tick_params(labelsize=axistick)

    ax0.set_ylabel('Depth (m)',fontsize=axisfont)
    ax1.set_ylabel('Depth (m)',fontsize=axisfont)
    ax2.set_ylabel('Depth (m)',fontsize=axisfont)
    ax3.set_ylabel('Depth (m)',fontsize=axisfont)

    ax1.set_xlabel('meters from outfall',fontsize=axisfont)
    ax3.set_xlabel('meters from outfall',fontsize=axisfont)

    ax0.tick_params(axis='both',which='major',labelsize=axistick)
    ax1.tick_params(axis='both',which='major',labelsize=axistick)
    ax2.tick_params(axis='both',which='major',labelsize=axistick)
    ax3.tick_params(axis='both',which='major',labelsize=axistick)

    #fig0.savefig(savename0+'%02d'%t_i+'.png',bbox_inches='tight')
    #fig2.savefig(savename2+'%02d'%t_i+'.png',bbox_inches='tight')
    fig0.savefig(savename0+'%02d'%(t_i+24)+'.png',bbox_inches='tight')
    fig2.savefig(savename2+'%02d'%(t_i+24)+'.png',bbox_inches='tight')
    plt.close('all')

