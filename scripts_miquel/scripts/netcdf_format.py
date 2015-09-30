import numpy as np
#import netCDF4 as nc4
from netCDF4 import Dataset
import os
import timeit
import datetime


# Originary netCDF files
fields_in = '/predator/emason/runs2009/gc_2009_1km_60/roms_avg.0990.nc'
grid_in = '/predator/emason/runs2009/gc_2009_1km_60/gc_2009_1km_grd_smooth.nc'

old_fields = Dataset(fields_in)
old_grid = Dataset(grid_in)

oldfiles = [old_fields, old_grid]

# Create the new netCDF files that tracpy will read
folder = 'nc_files_for_tracpy2'
if not os.path.exists(folder):
	os.makedirs(folder)
new_fields = Dataset(folder + '/roms_fileds.nc', 'w', format='NETCDF3_64BIT')
new_grid = Dataset(folder + '/roms_grid.nc', 'w', format='NETCDF3_64BIT')

newfiles = [new_fields, new_grid]

# Copy all attributes from originary to new files

'''new_fields.CPPS = old_fields.CPPS
new_fields.SRCS = old_fields.SRCS
new_fields.VertCoordType = old_fields.VertCoordType
#new_fields.cmptypes = old_fields.cmptypes
#new_fields.data_model = old_fields.data_model
#new_fields.dimensions = old_fields.dimensions
#new_fields.disk_format = old_fields.disk_format
new_fields.dt = old_fields.dt
new_fields.dtfast = old_fields.dtfast
new_fields.file_format = old_fields.file_format
new_fields.gamma2 = old_fields.gamma2
new_fields.grid_file = new_grid
new_fields.groups = old_fields.groups
new_fields.init_file = old_fields.init_file
new_fields.keepweakref = old_fields.keepweakref
new_fields.navg = old_fields.navg
new_fields.ndtfast = old_fields.ndtfast
new_fields.ntimes = old_fields.ntimes
new_fields.ntsavg = old_fields.ntsavg
new_fields.nwrt = old_fields.nwrt
new_fields.parent = old_fields.parent
new_fields.path = old_fields.path
new_fields.rdrg = old_fields.rdrg
new_fields.rdrg2 = old_fields.rdrg2
new_fields.rdrg2_units = old_fields.rdrg2_units
new_fields.rdrg_units = old_fields.rdrg_units
new_fields.rho0 = old_fields.rho0
new_fields.title = old_fields.title
new_fields.tnu2 = old_fields.tnu2
new_fields.tnu2_units = old_fields.tnu2_units
new_fields.type = old_fields.type
new_fields.variables = old_fields.variables
new_fields.visc2 = old_fields.visc2
new_fields.visc2_units = old_fields.visc2_units
new_fields.vltypes = old_fields.vltypes

new_grid.VertCoordType = old_grid.VertCoordType
#new_grid.cmptypes = old_grid.cmptypes
new_grid.data_model = old_grid.data_model
new_grid.date = old_grid.date
new_grid.dimensions = old_grid.dimensions
new_grid.disk_format = old_grid.disk_format
new_grid.file_format = old_grid.file_format
new_grid.groups = old_grid.groups
new_grid.keepweakref = old_grid.keepweakref
new_grid.parent = old_grid.parent
new_grid.path = old_grid.path
new_grid.title = old_grid.title
new_grid.type = old_grid.type
new_grid.variables = old_grid.variables
new_grid.vltypes = old_grid.vltypes'''

t = datetime.datetime(2015, 7, 24, 16, 0)
i = 0
for element in newfiles:
	print element, oldfiles[i]
	# Copy dimensions
	print '----Dimensions:'
	for dname, the_dim in oldfiles[i].dimensions.iteritems():
		print dname, len(the_dim), the_dim
		element.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
	# Copy variables
	print '----Variables:'
	for v_name, varin in oldfiles[i].variables.iteritems():
		outVar = element.createVariable(v_name, varin.datatype, varin.dimensions)
		print v_name, varin.datatype, varin.dimensions

		# Copy variable attributes
		outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
		print('    ','Variable attributes copied. Number of dictionaries copied: ',
					len(outVar.ncattrs()))
		print '	', (datetime.datetime.now()-t).total_seconds()
			#print('	Time needed: ', timeit.timeit('outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})'))
		outVar[:] = varin[:]
		print('    ','Variable values copied. Nuumber of values copied: ', 
					outVar[:].size)
		print '	', (datetime.datetime.now()-t).total_seconds()
			#print('	Time needed: ', timeit.timeit('outVar[:] = varin[:]'))
	# Close output file
	#element.close()

	i = i + 1

print 'Out of the copy loop. Add some more stuff.'

# Now it is time to modify the new files to give them the format we want
print 'Now it is time to modify the new files to give them the format we want.'
print 'Adding variables to new_fields:'
# s_w
print 's_w'
s_w = old_fields.dimensions['s_w']
#new_fields.createDimension('s_w', len(s_w))
Nlev = len(s_w)
nf_s_w = new_fields.createVariable('s_w', 'd', 's_w')
nf_s_w[:] = np.arange(-1, 1. / Nlev, 1. / (Nlev-1), dtype='d')[:]
'''for i in range(len(s_w)):
	nf_s_w[i] = i'''

# Cs_w
print 'Cs_w'
nf_Cs_w = new_fields.createVariable('Cs_w', 'float32', 's_w')
nf_Cs_w[:] = old_fields.Cs_w[:]

# hc
print 'hc'
nf_hc = new_fields.createVariable('hc', 'float32', dimensions=())
nf_hc[:] = old_fields.hc

# theta_s
print 'theta_s'
nf_theta_s = new_fields.createVariable('theta_s', 'float32', dimensions=())
nf_theta_s[:] = old_fields.theta_s

# theta_b
print 'theta_b'
nf_theta_b = new_fields.createVariable('theta_b', 'float32', dimensions=())
nf_theta_b[:] = old_fields.theta_b

print 'Adding variables to new_grids:'
# lon and lat u and v
print 'lon and lat u and v'
lon_rho = old_grid.variables['lon_rho'][:]
lat_rho = old_grid.variables['lat_rho'][:]

xi_u = old_fields.dimensions['xi_u']
eta_v = old_fields.dimensions['eta_v']
new_grid.createDimension('xi_u', len(xi_u))
new_grid.createDimension('eta_v', len(eta_v))

'''xi_rho = old_fields.dimensions['xi_rho']
eta_rho = old_fields.dimensions['eta_rho']
new_grid.createDimension('xi_rho', len(xi_rho))
new_grid.createDimension('eta_rho', len(eta_rho))'''

ng_lon_u = new_grid.createVariable('lon_u', 'float64', ('eta_v', 'xi_rho'))
ng_lat_u = new_grid.createVariable('lat_u', 'float64', ('eta_v', 'xi_rho'))
ng_lon_v = new_grid.createVariable('lon_v', 'float64', ('eta_rho', 'xi_u'))
ng_lat_v = new_grid.createVariable('lat_v', 'float64', ('eta_rho', 'xi_u'))

Mp, Lp = lon_rho.shape
M = Mp-1
L = Lp-1

ng_lon_u[:, :] = 0.5*(lon_rho[0:M, :] + lon_rho[1:Mp, :])
ng_lat_u[:, :] = 0.5*(lat_rho[0:M, :] + lat_rho[1:Mp, :])
ng_lon_v[:, :] = 0.5*(lon_rho[:, 0:L] + lon_rho[:, 1:Lp])
ng_lat_v[:, :] = 0.5*(lat_rho[:, 0:L] + lat_rho[:, 1:Lp])


new_fields.close()
new_grid.close()

print 'Done.'
