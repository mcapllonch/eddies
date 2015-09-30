from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors, cbook
from matplotlib.colors import LightSource
import sys
import os
from math import exp
from datetime import datetime
import gsw
import numexpr as ne
import octant
from scipy.interpolate import griddata
from PIL import Image


def vorticity(u, v, dx, dy):
    '''
    '''
    if u.ndim > 2: # 3d
        _, uy, _ = np.gradient(u, 0, dx, dy)
        _, _, vx = np.gradient(v, 0, dx, dy)
    else: # 2d
        uy, _ = np.gradient(u, dx, dy)
        _, vx = np.gradient(v, dx, dy)

    return ne.evaluate('vx - uy')


def interp_field_1D(field, nnew, axis=0, which=None):
	## 'which' option allows to create an array nnew times smaller, thus using a
	## lot less memory.
	## Create the new array for the field
	ns = []
	for i, val in enumerate(field.shape):
		if i == axis:
			if which is None:
				ns.append(val*nnew)
			else:
				ns.append(1)
		else:
			ns.append(val)		
	newshape = tuple(ns)
	newf = np.zeros(newshape)

	## Start interpolation:
	#for ax in range(field.ndim):
	newfsh0 = newf.shape[0]
	if axis == 0:
		try:
			for oldf_ti in range(field.shape[0]-1):
				df = field[oldf_ti + 1, :] - field[oldf_ti, :]
				newf[:] = field[oldf_ti, :] + (which/nnew)*df[:]
			#return newf
		except NameError:
			for oldf_ti in range(field.shape[0]-1):
				df = field[oldf_ti + 1, :] - field[oldf_ti, :]
				for n in range(nnew):
					newf[oldf_ti*nnew + n, :] = field[oldf_ti, :] + (n/nnew)*df[:]
		return newf
	else:
		print 'axis =! 0 not implemented yet. Leaving the field as is.'
		return


ls = LightSource(0, -60)		
bath_res_fac = 2
nfcontours = 60
lw = 8.

cmap1 = cm.gist_stern_r
cmap2 = cm.bone
cmap3 = cm.seismic

# Camera
move_cam = False
height = 40
azimut = -140


t0 = 199
tsep = 1
dpi = 1000
anim_nu = 140 # If it is 139, start from frame 14.
anim_name = 'Animation_%04i'%anim_nu
folder_path = '/home/mcapllonch/Documents/Animacions/' + anim_name + '/'

if not os.path.exists(folder_path):
	os.makedirs(folder_path)

roms_file = '/home/mcapllonch/scripts_miquel/nc_files_for_tracpy/roms_fileds.nc'
	# Grids
grd_roms3 = '/home/mcapllonch/scripts_miquel/nc_files_for_tracpy/roms_grid.nc'

# Read grid variables and declare variables allocating memory for them
with Dataset(grd_roms3) as nc:
	x = nc.variables['lon_rho'][:]
	y = nc.variables['lat_rho'][:]
	h_neg = -nc.variables['h'][:]
	mask_rho = nc.variables['mask_rho'][:]
with Dataset(roms_file) as nc:
	roms_times = nc.variables['ocean_time'].size
	var_shape = nc.variables['temp'][0, -1].shape
	u_shape = nc.variables['u'][0, -1].shape
	v_shape = nc.variables['v'][0, -1].shape


mask_rho_2 = np.array([mask_rho, mask_rho])
pt = np.ma.array(np.zeros((2, ) + var_shape), mask=mask_rho_2)
zeta = np.ma.array(np.empty_like(pt), mask=mask_rho_2)
#u = np.ma.array(np.zeros((2, ) + u_shape), mask=mask_rho_2)
#v = np.ma.array(np.zeros((2, ) + v_shape), mask=mask_rho_2)
#u = np.zeros((2, ) + u_shape)
#v = np.zeros((2, ) + v_shape)
u = np.zeros((2, ) + var_shape)
v = np.zeros((2, ) + var_shape)

pt_int = np.zeros(var_shape)
zeta_int = np.empty_like(pt_int)

print 'Calculating means, maxima and minima'

pt_means = np.zeros(roms_times)
pt_maxims = np.zeros(roms_times)
pt_minims = np.zeros(roms_times)

zeta_means = np.zeros(roms_times)
zeta_maxims = np.zeros(roms_times)
zeta_minims = np.zeros(roms_times)

for t_r in np.arange(roms_times):
	print t_r, ' of ', roms_times	
	## Read Fields
	with Dataset(roms_file) as nc:
		zeta[0] = nc.variables['zeta'][t_r]
		pt[0] = nc.variables['temp'][t_r, -1]
		pt[np.where(pt == pt.min())] = pt[np.where(pt != pt.min())].min()

	pt_means[t_r] = pt[0].mean()
	pt_maxims[t_r] = pt[0].max()
	pt_minims[t_r] = pt[0].min()

	zeta_means[t_r] = zeta[0].mean()
	zeta_maxims[t_r] = zeta[0].max()
	zeta_minims[t_r] = zeta[0].min()

pt_mean = pt_means.mean()
pt_maxim = pt_maxims.max()
pt_minim = pt_minims.min()

zeta_mean = zeta_means.mean()
zeta_maxim = zeta_maxims.max()
zeta_minim = zeta_minims.min()

'''# Mask 0 Celsius values for temperature
combined_mask = np.ma.mask_or(mask_rho_2.astype('bool'), 
	(np.ma.masked_less(pt, pt[np.where(pt != pt.min())].min()).mask))
pt = np.ma.array(pt.data, mask=combined_mask)'''

a = np.linspace(24., pt_maxim, 40)
b = np.linspace(22.5, a.min(), 30)
c = np.linspace(pt_minim, b.min(), 10)
levels_pt = np.concatenate((c, b[1:], a[1:]))
#levels_pt = np.concatenate((np.linspace(pt_minim, a.min(), 10), a)).flatten()

#levels_pt = np.linspace(pt_minim + 0.1, pt_maxim, nfcontours)
#levels_pt = -np.cos(levels_pt * np.pi / pt_maxim) * pt_maxim

# Loop
t_roms = -1
frame = t0
nsteps = 10

ntimes = roms_times * nsteps
# Camera
height_cam = (np.arange(ntimes)/float(ntimes) + 0.5) * height
azimut_cam = (np.arange(ntimes)/float(ntimes)) * azimut

for time_p in np.arange(t0, ntimes, tsep):

	t1 = datetime.now()

	t_roms_prev = t_roms
	t_roms = time_p // nsteps # Current step for the ROMS fields.

	t4int = time_p # Current step for the interpolated ROMS fields.
	
	## Read Fields
	read_now = t_roms != t_roms_prev
	if read_now:
		print 'Reading fields from ROMS', t_roms
		with Dataset(roms_file) as nc:
			pt[0] = nc.variables['temp'][t_roms, -1]
			pt[1] = nc.variables['temp'][t_roms + 1, -1]
			u[0, :, :-1] = nc.variables['u'][t_roms, -1]
			v[0, :-1, :] = nc.variables['v'][t_roms, -1]	
			u[1, :, :-1] = nc.variables['u'][t_roms + 1, -1]
			v[1, :-1, :] = nc.variables['v'][t_roms + 1, -1]
			zeta[0] = nc.variables['zeta'][t_roms]
			zeta[1] = nc.variables['zeta'][t_roms + 1]

		vort = vorticity(u[0], v[0], 1.e3, 1.e3)
		#pt[np.where(pt == pt.min())] = pt[np.where(pt != pt.min())].min()
		mask_zeros = np.ma.masked_less(pt[0], pt_minim).mask
		combined_mask = np.ma.mask_or(~mask_rho_2.astype('bool')[0], mask_zeros)
		pt[0] = np.ma.array(pt[0].data, mask=combined_mask)
		pt[0].data[np.where(pt[0].mask == True)] = np.nan

		## Trampa:
		'''pt.data[np.where(pt.data <= 22.)] = 22.
		levels_pt = np.linspace(21.9, pt_maxim, nfcontours)'''

	## Interpolate fields in time

	pt_int[:] = interp_field_1D(pt[:], nsteps, axis=0, which=time_p-nsteps*t_roms)[:]
	zeta_int[:] = interp_field_1D(zeta[:], nsteps, axis=0, which=time_p-nsteps*t_roms)[:]

	#fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
	fig = plt.figure(figsize=(19.2, 10.8))
	ax = Axes3D(fig)
	m = Basemap(llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(),)

	#im2 = ax.plot_surface(x, y, h_neg, cmap=cm.coolwarm)
	## Plot
	rgb = ls.shade(h_neg[::bath_res_fac, ::bath_res_fac], cmap=cmap2)
	im2 = ax.plot_surface(x[::bath_res_fac, ::bath_res_fac], y[::bath_res_fac, ::bath_res_fac], 
							h_neg[::bath_res_fac, ::bath_res_fac], rstride=1, cstride=1, facecolors=rgb,
	                       linewidth=0., antialiased=False, shade=False,
	                       alpha=1.)


	cset2 = ax.contour(x, y, pt_int, levels=levels_pt, linewidth=lw, zdir='z', offset=0., cmap=cmap1, alpha=.8)
	#cset2 = ax.contourf(x, y, pt[0], levels=levels_pt, zdir='z', offset=0., cmap=cmap1, alpha=0.9, antialiased=True)
	#cset3 = ax.contour(x, y, h_neg, levels=np.array([-10000, 0, 10000]), zdir='z', offset=0., colors='k')
	#cset3 = ax.contour(x, y, -h_neg, levels=np.linspace(-h_neg.max() + 0.1, -h_neg.max(), 3), colors='k', linewidth=20.)
	#ax.add_collection3d(m.drawcoastlines(linewidth=0.25))
	cset4 = ax.contour(x, y, mask_rho, colors='k', linewidths=1.)
	cset5 = ax.contourf(x, y, mask_rho, levels=np.array([-0.5, 0., 0.5]), cmap=cm.binary)
	#cset5 = ax.contourf(x, y, mask_rho, levels=np.array([-0.5, 0., 0.5]), colors='#D6C6A9')
	#cset6 = ax.contour(x, y, zeta[0], levels=np.linspace(-0.3, 0.3, 20), colors='k')
	#cset6 = ax.contour(x, y, zeta_int, levels=np.linspace(-0.3, 0., 20), linewidths=2., colors='k')
	#cset7 = ax.contour(x, y, zeta_int, levels=np.linspace(0.18, 0.5, 20), linewidths=2., colors='k')
	#cset5 = ax.contour(x, y, vort, nfcontours, cmap=cmap3)
	#m.bluemarble()

	if move_cam:
		ax.view_init(height_cam[time_p], azimut_cam[time_p])
	else:
		ax.view_init(height, azimut)
	plt.axis('off')

	image_name = folder_path + anim_name + '_Frame_%04i'%frame + '.png'
	fig.savefig(image_name, dpi=dpi)
	#plt.show()
	plt.close(fig)
	## Resizing a large image to the wanted resolution
	size = 1920, 1080
	image = Image.open(image_name)
	im_resized = image.resize(size, Image.ANTIALIAS)
	im_resized.save(image_name, "PNG")

	t2 = datetime.now()
	print 'Frame created.', anim_nu, frame, (t2 - t1).total_seconds()
	frame += 1