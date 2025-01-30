from netCDF4 import Dataset
import numpy as np

nc = Dataset('fb_edata.nc','w')

nc.createDimension('two',2)
unc = nc.createVariable('u',np.double,('two'))
vnc = nc.createVariable('w',np.double,('two'))

unc.output_period = '10'
