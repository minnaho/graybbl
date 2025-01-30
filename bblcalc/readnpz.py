import numpy as np
import matplotlib.pyplot as plt

plt.ion()



# And I just realized I gave you the wrong PGF for the no rotation case, that one is:
#PGF=2x10^-7
datastr = dict(np.load('flat_BBL_f_strat.npz',allow_pickle=True))
# PGF = u_g*f = 6x10^-6
datawif = dict(np.load('flat_BBL_f.npz',allow_pickle=True))

wif_uprof = datawif['u_prof']
wif_corfq = datawif['cor_freq']
wif_dpdx  = datawif['dpdx_prof']
wif_avgdpdx = datawif['avg_dpdx']
wif_dz    = datawif['dz']
wif_zplt  = wif_dz*np.arange(1,datawif['nz']+wif_dz)


str_dz    = datastr['dz']
str_dx    = datastr['dx']
str_zplt  = str_dz*np.arange(1,datastr['nz']+str_dz)
str_xplt  = str_dx*np.arange(1,datastr['nx']+str_dx)
str_tplt  = np.arange(1,datastr['u_prof'].shape[0]+1)*5 # time in minutes
str_avg_u = datastr['avg_u']
str_uprof = datastr['u_prof']

fig1,ax1 = plt.subplots(1,3,figsize=[15,6])
ax1.flat[0].plot(wif_uprof,wif_zplt)
ax1.flat[0].set_xlabel('u')
ax1.flat[1].plot(wif_dpdx,wif_zplt)
ax1.flat[1].set_xlabel('dpdx')
ax1.flat[2].plot(wif_avgdpdx,wif_zplt)
ax1.flat[2].set_xlabel('avg_dpdx')

# LES
fig2,ax2 = plt.subplots(2,1,figsize=[15,10])
# last time step of u (averged over last 5 minutes)
avguplt2 = ax2.flat[0].pcolormesh(str_xplt,str_zplt,str_avg_u[-1],cmap='rainbow',vmin=0,vmax=0.061)
cb2ax0 = fig2.colorbar(avguplt2,ax=ax2.flat[0])
uprofplt2 = ax2.flat[1].pcolormesh(str_tplt,str_zplt,str_uprof.transpose(),cmap='rainbow',vmin=0,vmax=0.061)
cb2ax1 = fig2.colorbar(uprofplt2,ax=ax2.flat[1])




