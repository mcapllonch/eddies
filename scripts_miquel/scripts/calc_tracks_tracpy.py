import numpy as np
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from tracpy.tracpy_class import Tracpy
from tracpy.tracpy_class import datetime

import matplotlib.pyplot as plt 
from matplotlib.mlab import find

# nc files
loc = '/home/mcapllonch/scripts_miquel/nc_files_for_tracpy/roms_fileds.nc'
grid_filename = '/home/mcapllonch/scripts_miquel/nc_files_for_tracpy/roms_grid.nc'

ndays = 29
#date = datetime.datetime(2009, 11, 25, 0)
date = datetime.datetime(2006, 3, 26, 0, 14, 56)
tseas = 24*3600
time_units = 'seconds since 1970-01-01'
nsteps = 48
N = 48
ff = 1
ah = 0. # m^2/s
av = 0. # m^2/s
doturb = 0
name = 'prova_24'
# Input starting locations as real space lon,lat locations
lon0, lat0 = np.meshgrid(np.linspace(-16.5, -16., 100), \
                            np.linspace(26.5, 27., 100))
#lon0 = np.linspace(-17., -15., 2)
#lat0 = 

#print 'lon i lat originals:', lon0, lat0
do3d = 1
#z0 = 's'
#z0 = np.zeros(lon0.size)
z0 = np.random.rand(lon0.size)
#z0 = np.ones(lon0.size)
#z0 = np.linspace(0, 1, lon0.size)
z0[:] = -120. + 10.*z0[:]
#exit()
#lon0orig = lon0
'''num_layers = 30
zpar = num_layers-1'''
zpar='fromZeta'
print 'ok 1'
tp = Tracpy(currents_filename=loc, grid_filename=grid_filename, name=name, 
			tseas=tseas, ndays=ndays, nsteps=nsteps, usebasemap=True, N=N, ff=ff,
			 ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, 
			 time_units=time_units)
print 'ok 2'
netCDF._netCDF4._set_default_format(format='NETCDF3_64BIT')
tp._readgrid()
print 'ok 3'

## Remove values of z0 that are deeper than the deepest level
from scipy.interpolate import griddata
import octant

h = tp.grid['h'].T.copy(order='c')
# Use octant to calculate depths for the appropriate vertical grid parameters
# have to transform a few back to ROMS coordinates and python ordering for this
zwt = octant.depths.get_zw(tp.grid['Vtransform'], tp.grid['Vstretching'], tp.grid['km']+1, tp.grid['theta_s'], tp.grid['theta_b'], 
                h, tp.grid['hc'], zeta=0, Hscale=3)

H = zwt[0, :, :].T
#Hant = tp.grid['zwt0'][:,:,0]
#H = tp.zwt[:, :, 1, 0]
coord = np.array([tp.grid['lonr'].flatten(),tp.grid['latr'].flatten()])
coord = np.transpose(coord)
coord0 = np.transpose(np.array([lon0.flatten(),lat0.flatten()]))
toto = griddata(coord,H.flatten(),coord0)
step = 100 # This is the maximum depth diference between two nearby gridpoints that this method can stand.
b = np.where(tp.z0-toto > step)
z0_old = tp.z0
tp.z0 = tp.z0[b]
#print 'lon i lat aplanats:', lon0.flatten(), lat0.flatten()
lon0_chosen = lon0.flatten()[b]
lat0_chosen = lat0.flatten()[b]
lon0 = lon0_chosen
lat0 = lat0_chosen
'''print 'lon i lat escollits:', lon0_chosen, lat0_chosen
print 'z0: ', z0
print 'toto: ', toto
print 'tp.z0: ', tp.z0'''

'''for i in range(H.shape[0]):
	for j in range(H.shape[1]):
		ind = zwt[i, j, :] <= tp.z0.min()'''

# Eliminate points that are outside domain or in masked areas
lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)
# Resize tp.z0 or lon0 and lat0, depending on which is smaller
if tp.z0.size >= lon0.size:
	print 'tp.z0 reduced', tp.z0.size, lon0.size
	tp.z0 = tp.z0[0:lon0.size+1]
'''else:
	print 'lon and lat reduced', tp.z0.size, lon0.size
	lon0 = lon0[0:tp.z0.size]
	lat0 = lat0[0:tp.z0.size]'''


lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)
tracpy.plotting.tracks(lonp, latp, tp.name, tp.grid)