import os
import sys
sys.path.append('/data/project3/minnaho/global/')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset,num2date,date2num
import glob
from scipy.signal import savgol_filter
import ROMS_depths as depths
plt.ion()

datanc = Dataset('pipe_3m_dt1_nh_F1_noplume_his.20000101120000.nc','r')
gridnc = Dataset('pipe_3m_dt1_nh_F1_noplume_grd.nc','r')

zr = depths.get_zr_zeta_tind(datanc,gridnc,-1)

unc = np.squeeze(datanc.variables['u'])

cent = int(unc.shape[1]/2.)
nearend = int(unc.shape[2]-20)

uprof = unc[:,cent,nearend]

fig,ax = plt.subplots(1,1)
ax.plot(uprof,zr[:,cent,nearend])

