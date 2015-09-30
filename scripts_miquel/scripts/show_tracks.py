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


print
print '--------------------------------------'
print '### Program for showing the floats ###'
print '--------------------------------------'
print

## Functions


def include_bathymetry(x, y):
	
	print 'Including bathymetry.'

	t1 = datetime.now()


	surf_cmap = cm.gist_earth
	lev = np.linspace(h_neg_basin.min(), h_neg_basin.max(), 1000)
	norml = colors.BoundaryNorm(lev, 256)

	ls = LightSource(0, -75)
	rgb = ls.shade(h_neg[::bath_res_fac, ::bath_res_fac], cmap=surf_cmap)
	## Plot
	im2 = ax.plot_surface(x[::bath_res_fac, ::bath_res_fac], y[::bath_res_fac, ::bath_res_fac],
	                       h_neg[::bath_res_fac, ::bath_res_fac], rstride=1, cstride=1, facecolors=rgb,
	                       linewidth=0, antialiased=True, shade=False,
	                       alpha=1.)
	
	'''h_neg_basin[np.where(np.ma.getmask(h_neg_basin)==True)] = np.nan
	im2 = ax.plot_surface(x, y, h_neg_basin, rstride=1, cstride=1, norm=norml,
	                       linewidth=0, antialiased=False, shade=False,
	                       alpha=0.8)'''
	t2 = datetime.now()
	print 'Bathymetry finished. Time elapsed:', (t2 - t1).total_seconds()


def fade(x, y, lons, lats, depths, fvalue, fwidth, fade_field=None, fade_method='Masking'):

	sigma = fwidth
	fotf = np.zeros(lons.shape)
	dnearest = np.zeros(lons.shape)

	for i in range(lons.size):
		ijdnearest = nearest_2D_single(lons[i], lats[i], x, y)
		fotf[i] = fieldonthefloat(ijdnearest[:-1], depths[i], field=fade_field)
		dnearest[i] = ijdnearest[-1]


	'''ijnearest = nearest_2D(lons, lats, x, y)
	for i in range(lons.size):
		fotf[i] = fieldonthefloat(ijnearest, depths[i], field=fade_field[0, :])'''

	if fade_method == 'Alphas':
		alphas = np.zeros(lons.shape)
		for i in range(lons.size):
			alphas[i] = exp(-((fotf[i]-fvalue)/sigma)**2)
	elif fade_method == 'Masking':
		flims = (fvalue - fwidth, fvalue + fwidth)
		dist_mask = np.ma.masked_greater(dnearest, dmax).mask
		scatt_mask = np.ma.masked_outside(fotf, flims[0], flims[1]).mask
		combined_mask_5 = np.ma.mask_or(dist_mask, scatt_mask)
		lons = np.ma.array(lons, mask=combined_mask_5)
		#lons = np.ma.array(lons, mask=scatt_mask)
		alphas = np.ones(lons.shape)
	elif fade_method == 'Both':
		flims = (fvalue - fwidth, fvalue + fwidth)
		scatt_mask = np.ma.masked_outside(fotf, flims[0], flims[1]).mask
		#combined_mask_5 = np.ma.mask_or(lons.mask, scatt_mask)
		#lons = np.ma.array(lons.data, mask=combined_mask_5)
		lons = np.ma.array(lons, mask=scatt_mask)

		alphas = np.zeros(lons.shape)
		for i in range(lons.size):
			alphas[i] = exp(-((fotf[i]-fvalue)/sigma)**2)

		print 'good lons:', lons.count()
		print alphas.min(), alphas.max()
		'''	elif: fade_method == 'Coloring':
		'''
	else:
		alphas = np.ones(lons.shape)

	if fade_outs_custom_reg:
		dist_mask = np.ma.masked_greater(dnearest, dmax).mask
		lons = np.ma.array(lons, mask=dist_mask)
		#alphas = np.ones(lons.shape)

	return alphas, lons


def trackscatter_all(lons, lats, depths, x, y, firstbar_fl, colorby=None, color='Green', fade_field=None, 
	colorfield=None, fade_method=None, what2do='scatter', zorder=None):

	## Select transparency if fade
	tt1 = datetime.now()
	if fade_method != None:
		alphas, lons = fade(x, y, lons, lats, depths, fvalue = fvalue, fwidth = fwidth, fade_field=fade_field, fade_method=fade_method)
	else:
		alphas = np.ones(lons.shape)

	tt2 = datetime.now()

	#print 'TIME for the fading process: ', (tt2 - tt1).total_seconds()

	## Select coloring method and plot
	colorplot = np.zeros(lons.size)
	alpha = 1 # To be modified if there is fading
	if colorby == 'depth':
		colorplot = depths
		im = ax.scatter(lons, lats, depths, marker='o', s=parts_size**2, c=colorplot, 
			cmap=cm.Purples, alpha=alpha, zorder=zorder)
		if firstbar_fl:
			cbar = fig.colorbar(im, orientation=cb_orientation, shrink=cb_shrink)
			cbar.set_label('Depth bar. Units = Unknown')
	elif colorby == 'day':
		im = ax.scatter(lons, lats, depths, c=colorplot, cmap=cm.Purples, alpha=alpha, zorder=zorder)
		if firstbar_fl:
			cbar = fig.colorbar(im, orientation=cb_orientation, shrink=cb_shrink)
			cbar.set_label('Time bar. Units = Days')
	elif colorby == 'field':
		for i in range(lons.size):
			ijdnearest = nearest_2D_single(lons[i], lats[i], x, y)
			colorplot[i] = fieldonthefloat(ijdnearest[:-1], depths[i], field=fade_field[0, :])
		#colorplot = colorplot/colorplot.max()
		im = ax.scatter(lons, lats, depths, c=colorplot, cmap=cm.jet, alpha=alpha, zorder=zorder)
		if firstbar_fl:
			cbar = fig.colorbar(im, orientation=cb_orientation, shrink=cb_shrink)
			cbar.set_label(labels[fill_field])
	elif colorby == 'fixed':

		im = ax.scatter(lons, lats, depths, c=color, alpha=alpha, zorder=zorder)

	else:
		#print 'In trackscatter_all. I am here.'
		if what2do == 'plot':
			tt3 = datetime.now()
			im = ax.plot(lons, lats, depths, marker='o', markersize=parts_size, lw=0, color=parts_color, markeredgecolor=parts_edge_col, 
				markeredgewidth=patrs_edge_size)
		elif what2do == 'scatter':
			tt3 = datetime.now()
			'''rgba_colors = np.zeros((lons.size, 4))
			rgba_colors[:, 3] = alphas
			#im = ax.scatter(lons, lats, depths, marker='o', s=2.5, color=rgba_colors)
			im = ax.scatter(lons, lats, depths, marker='o', s=1.5, color=rgba_colors)'''
			im = ax.scatter(lons, lats, depths, marker='o', s=parts_size**2, color=parts_color, edgecolors=parts_edge_col, linewidth=patrs_edge_size)


	firstbar_fl = False


	'''def show_slice_imit_isosurfs(x, y, var, var_raw, level, firstbar_fie, cmap, alpha, label, fill=False, bar=False, 
	colors=None, pos_neg=False, zorder=None):
	if pos_neg:
		abv = abs(var_raw)
		lims = (-abv[np.where(~(var_raw == var_raw.min()))].max(), abv[np.where(~(var_raw == var_raw.min()))].max())
		#lims = (- 15., 15.)
		#lims = (fvalue - fwidth, fvalue + fwidth)
	else:
		lims = (var_raw[np.where(~(var_raw == var_raw.min()))].min(), var_raw[np.where(~(var_raw == var_raw.min()))].max())
	if fill:
		cset = ax.contourf(x, y, var, zdir='z', offset=level, 
			levels=np.linspace(lims[0], lims[1], nfcontours), cmap=cmap, 
			alpha=alpha, antialiased=True, zorder=zorder)
	else:
		if colors == None:
			cset = ax.contour(x, y, var, zdir='z', offset=level, 
				#levels=np.linspace(var_raw.min(), var_raw.max(), 100), cmap=cm.jet, 
				levels=np.linspace(lims[0], lims[1], ncontours), cmap=cmap, 
				alpha=alpha, zorder=zorder, linewidth=cont_lw)
		else:
			cset = ax.contour(x, y, var, zdir='z', offset=level, 
				#levels=np.linspace(var_raw.min(), var_raw.max(), 100), cmap=cm.jet, 
				levels=np.linspace(lims[0], lims[1], ncontours), colors=colors, 
				alpha=alpha, zorder=zorder, linewidth=cont_lw)

	if firstbar_fie and bar:
	#if bar:
		cbar = fig.colorbar(cset, orientation=cb_orientation, shrink=cb_shrink) #, cax=fig.add_axes([0.8, 0.1, 0.03, 0.8]))
		cbar.set_label(label)
		ticks = np.linspace(lims[0], lims[1], 11)
		cbar.ax.set_yticks(ticks)
		cbar.ax.set_yticklabels(['%0.1f'%tick for tick in ticks])

	firstbar_fie = False'''


def show_slice_imit_isosurfs2(x, y, var, vmin, vmax, level, firstbar_fie, cmap, alpha, label, fill=False, bar=False, 
	colors=None, pos_neg=False, zorder=None):
	if pos_neg:
		abv[:] = np.abs(np.array([vmin, vmax]))
		lims = (-abv.max(), abv.max())
	else:
		lims = (vmin, vmax)
	if fill:
		cset = ax.contourf(x, y, var, zdir='z', offset=level, 
			levels=np.linspace(lims[0], lims[1], nfcontours), cmap=cmap, 
			alpha=alpha, antialiased=True, zorder=zorder)
	else:
		if colors == None:
			cset = ax.contour(x, y, var, zdir='z', offset=level, 
				levels=np.linspace(lims[0], lims[1], ncontours), cmap=cmap, 
				alpha=alpha, zorder=zorder, linewidth=cont_lw)
		else:
			cset = ax.contour(x, y, var, zdir='z', offset=level, 
				levels=np.linspace(lims[0], lims[1], ncontours), colors=colors, 
				alpha=alpha, zorder=zorder, linewidth=cont_lw)

	if firstbar_fie and bar:
	#if bar:
		cbar = fig.colorbar(cset, orientation=cb_orientation, shrink=cb_shrink) #, cax=fig.add_axes([0.8, 0.1, 0.03, 0.8]))
		cbar.set_label(label)
		ticks = np.linspace(lims[0], lims[1], 11)
		cbar.ax.set_yticks(ticks)
		cbar.ax.set_yticklabels(['%0.1f'%tick for tick in ticks])

	firstbar_fie = False


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


def nearest_2D_single(xpoint, ypoint, x, y):
	## Find the closest point in the grid to any given point.
	## Only in 2D. Not implemented in the vertical.

	d = np.hypot(x - xpoint, y - ypoint)
	inearest, jnearest = np.unravel_index(d.argmin(), d.shape)

	try:
		_ = inearest
	except:
		print 'nearest_2D failure.', i, j, xpoint, ypoint, x[i, j], y[i, j], d, dold
	return inearest, jnearest, d[inearest, jnearest]


def surrounding4(xpoint, ypoint, x, y):
	"""
	Find the 4 surrounding points in the grid to any given point.
	Only in 2D. Not implemented in the vertical.

	xpoint and ypoint: Single lonlat point.

	WARNING: I realized that this function does not really get the 4 surrounding oints, 
	but the 4 closest points (when I wrote it I thought they were going to be the same, 
	but they are not in this case, given the latlon distribution in the matrix), so it 
	would be better to use the function made by Evan Mason instead.
	"""

	d = np.hypot(x - xpoint, y - ypoint)

	i, j = (np.zeros(4, dtype=int), np.zeros(4, dtype=int))
	for ii in np.arange(4):
		i[ii], j[ii] = np.unravel_index(d.argmin(), d.shape)
		d[i[ii], j[ii]] = d.max()

	try:
		_ = i
	except:
		print 'surrounding4 failure.', i, j, xpoint, ypoint, x[i, j], y[i, j], d, dold
	return i, j


def nearest_2D(lon_pts, lat_pts, x, y):
	## Find the closest point in the grid to any given point.
	## Only in 2D. Not implemented in the vertical.

	ls = lon_pts.size
	xx = np.tile(x, (lon_pts.size, 1)).reshape((ls, x.shape[0], x.shape[1]))
	yy = np.tile(y, (lon_pts.size, 1)).reshape((ls, x.shape[0], x.shape[1]))
	lon_pts = np.repeat(lon_pts, x.size).reshape(xx.shape)
	lat_pts = np.repeat(lat_pts, x.size).reshape(xx.shape)

	d = np.hypot(xx - lon_pts, yy - lat_pts)

	ij = np.zeros(ls)
	for i in range(ls):
		ij[i] = np.unravel_index(d[i, :].argmin(), d.shape[1:])

	try:
		_ = inearest
	except:
		print 'nearest_2D failure.', i, j, xpoint, ypoint, x[i, j], y[i, j], d, dold
	return ij


def fieldonthefloat(ij, depth, field):
	levlist = list(lev_fie_chosen)
	levlist[0] = 10
	for l, prof in enumerate(levlist[:-1]):
		if depth == l:
			fint = field[l, ij[0], ij[1]]
			#print 'depth ok for l =', l
		elif (depth < prof) and (depth > levlist[l+1]):
			#print field.shape, l, ij[0], ij[1]
			fm1 = field[l, ij[0], ij[1]]
			fp1 = field[l+1, ij[0], ij[1]]
			d = prof - levlist[l+1]
			fint = fm1 + ((prof - depth)/d)*(fp1 - fm1)
			#print 'depth ok for l =', l
		#else:
			#print l, prof, depth, levlist[l+1]
	try:
		_ = fint
	except:
		print 'fieldonthefloat failure.', levlist, depth, l, prof
	return fint


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


def rho_eos_duko(temp, salt, z_r, rho0=1027.40898956694):
    '''
    Get Dukowicz density
    Input: rho_eos_duko(temp, salt, z_r, rho0=1027.40898956694)
    Output: rho
    '''
    temp = np.float64(temp)
    salt = np.float64(salt)
    rho0 = rho0
    temp0 = temp0
    salt0 = salt0

    r00 = 999.842594;   r01 = 6.793952E-2;  r02 = -9.095290E-3
    r03 = 1.001685E-4;  r04 = -1.120083E-6; r05 = 6.536332E-9

    r10 = 0.824493;     r11 = -4.08990E-3;  r12 = 7.64380E-5
    r13 = -8.24670E-7;  r14 = 5.38750E-9

    rS0 = -5.72466E-3;  rS1 = 1.02270E-4;   rS2 = -1.65460E-6
    r20 = 4.8314E-4

    K00 = 19092.56;     K01 = 209.8925;     K02 = -3.041638
    K03 = -1.852732e-3; K04 = -1.361629e-5

    K10 = 104.4077;     K11 = -6.500517;    K12 = 0.1553190
    K13 = 2.326469e-4

    KS0 = -5.587545;    KS1 = +0.7390729;   KS2 = -1.909078e-2

    #B00 = 0.4721788;    B01 = 0.01028859;   B02 = -2.512549e-4
    #B03 = -5.939910e-7
    #B10 = -0.01571896;  B11 = -2.598241e-4; B12 = 7.267926e-6
    #BS1 = 2.042967e-3

    #E00 = 1.045941e-5;  E01 = -5.782165e-10; E02 = 1.296821e-7
    #E10 = -2.595994e-7; E11 = -1.248266e-9;  E12 = -3.508914e-9

    sqrt_salt0 = np.sqrt(salt0)
    dr00 = r00 - 1000.e0
    #rho1_0 = dr00 + temp0 * (r01 + temp0 *(r02 + temp0 * (r03 + temp0 * (r04 + temp0 * r05)))) \
                            #+salt0*(r10+temp0*(r11+temp0*(r12+temp0*(    \
                                               #r13+temp0*r14)))          \
          #+sqrt_salt0*(rS0+temp0*(rS1+temp0*rS2))+salt0*r20)

    K0_Duk = temp0*( K01+temp0*( K02+temp0*( K03+temp0*K04 )))  \
          +salt0*( K10+temp0*( K11+temp0*( K12+temp0*K13 )) \
               +sqrt_salt0*( KS0+temp0*( KS1+temp0*KS2 )))

    dr00   = r00 - rho0
    sqrt_salt = ne.evaluate('sqrt(salt)')

    rho1 = ne.evaluate('( dr00+temp*( r01+temp*( r02+temp*( r03+temp*(          \
                                       r04+temp*r05 ))))       \
                      +salt*( r10+temp*( r11+temp*( r12+temp*( \
                                        r13+temp*r14 )))       \
                         +sqrt_salt*(rS0+temp*(                \
                               rS1+temp*rS2 ))+salt*r20 ))')

    K0 = ne.evaluate('temp * (K01 + temp * (K02 + temp*(K03 + temp * K04)))          \
      + salt * (K10 + temp * (K11 + temp * (K12 + temp * K13))           \
      + sqrt_salt * (KS0 + temp * (KS1 + temp * KS2)))')

    qp1 = ne.evaluate('0.1 * (rho0 + rho1) * (K0_Duk - K0) / ((K00 + K0) * (K00 + K0_Duk))')
    qp2 = 0.0000172
    # rho is returned by make_clim_ssh
    rho = ne.evaluate('rho1 + qp1 * abs(z_r) * (1. - qp2 * abs(z_r))')
    return rho, rho1, qp1, rho0 # changed 19/5/2014; before it returned rho + rho0


def omega2w(u_in, v_in, omega_in, zr_in, pm_in, pn_in):
    '''
    Compute vertical velocity ROMS output u,v,omega
   
        w = omega2w(u, v, omega, zr, pm, pn)

        Function by IMEDEA.
   
    '''
    import numexpr as ne
    
    u = u_in.view().swapaxes(2,0)
    v = v_in.view().swapaxes(2,0)
    omega = omega_in.view().swapaxes(2,0)
    zr = zr_in.view().swapaxes(2,0)
    pm = pm_in.view().swapaxes(1,0)
    pn = pn_in.view().swapaxes(1,0)
   
    Lr, Mr, Np = omega.shape
    N = Np - 1
    w   = np.zeros((Lr,   Mr,   N))
    Wxi = np.zeros((Lr-1, Mr,   N))
    Wet = np.zeros((Lr,   Mr-1, N))

    #  dsigma/dt component
    o_N = omega[:,:,N].view()
    o_Nm1 = omega[:,:,N-1].view()
    o_Nm2 = omega[:,:,N-2].view()
    w[:,:,N-1] = ne.evaluate('0.375 * o_N + 0.75 * o_Nm1 - 0.125 * o_Nm2')
    o_2N = omega[:,:,2:N].view()
    o_1Nm1 = omega[:,:,1:N-1].view()
    o_3Np1 = omega[:,:,3:N+1].view()
    o_Nm2 = omega[:,:,:N-2].view()
    w[:,:,1:N-1] = ne.evaluate('0.5625 * (o_2N + o_1Nm1) - 0.0625 * (o_3Np1 + o_Nm2)')
    o_2 = omega[:,:,2].view()
    o_1 = omega[:,:,1].view()
    o_0 = omega[:,:,0].view()
    w[:,:,0] = ne.evaluate('-0.125 * o_2 + 0.75 * o_1 + 0.375 * o_0')

    # domega/dx, domega/dy component
    pm_1L = pm[1:Lr].view()
    pm_Lm1 = pm[0:Lr-1].view()
    pn_1M = pn[:,1:Mr].view()
    pn_Mm1 = pn[:,0:Mr-1].view()
    for kk in np.nditer(np.arange(N)):
       
        u_kk = u[:,:,kk].view()
        zr_1L_kk = zr[1:Lr,:,kk].view()
        zr_Lm1_kk = zr[0:Lr-1,:,kk].view()
        #print 'u_kk, zr_1L_kk, zr_Lm1_kk, pm_1L, pm_Lm1: ', u_kk.shape, zr_1L_kk.shape, zr_Lm1_kk.shape, pm_1L.shape, pm_Lm1.shape
        Wxi[:Lr-1,:,kk] = ne.evaluate('(u_kk * (zr_1L_kk - zr_Lm1_kk)) * 0.5 * (pm_1L + pm_Lm1)')
       
        v_kk = v[:,:,kk].view()
        zr_1M_kk = zr[:,1:Mr,kk].view()
        zr_Mm1_kk = zr[:,0:Mr-1,kk].view()
        Wet[:,:Mr-1,kk] = ne.evaluate('(v_kk * (zr_1M_kk - zr_Mm1_kk)) * 0.5 * (pn_1M + pn_Mm1)')
   
    w_LM = w[1:Lr-1,1:Mr-1].view()
    Wxi_1 = Wxi[0:Lr-2,1:Mr-1].view()
    Wxi_2 = Wxi[1:Lr-1,1:Mr-1].view()
    Wet_1 = Wet[1:Lr-1,0:Mr-2].view()
    Wet_2 = Wet[1:Lr-1,1:Mr-1].view()
    w[1:Lr-1,1:Mr-1] =  ne.evaluate('w_LM + 0.5 * (Wxi_1 + Wxi_2 + Wet_1 + Wet_2)')

    w[:,0]    = w[:,1]
    w[:,Mr-1] = w[:,Mr-2]
    w[0]      = w[1]
    w[Lr-1]   = w[Lr-2]

    return w.view().swapaxes(2, 0)


def no_excursion(framesperday, maxdd):
	ntimes_p = lons_nc.shape[0]
	nparts = lons_nc.shape[1]
	convfact = 6.*framesperday/(maxdd)

	depths_nex = np.zeros(depths_nc.shape)
	depths_nex[:, :] = depths_nc[:, :]

	ddabs = np.zeros(nparts)
	ddnew = np.zeros(nparts)

	for time in range(1, ntimes_p - 1):
		inds = np.where(depths_nc.mask[time, :] == False)
		#inds2 = np.where(depths_nc.mask[time + 1, :] == False)
		ddabs[inds] = depths_nc[time + 1, inds] - depths_nc[time, inds] 
		ddnew[inds] = (1. - 0.99*(2./np.pi)*np.arctan(abs(ddabs[inds])*convfact))*ddabs[inds]
		depths_nex[time + 1, inds] = depths_nex[time, inds] + ddnew[inds]

		if time == 833:
			print 'In no_excursion.', ddabs[76], ddnew[76], depths_nc[time, 76], depths_nc[time + 1, 76], depths_nex[time, 76], depths_nex[time + 1, 76]
	depths_nex = np.ma.array(depths_nex, mask=depths_nc.mask)
	return depths_nex


def omega_to_w(omega, u_s, v_s, zeta, h_neg_, bound):

	# IMEDEA modification: omega is perpendicular to sigma surfaces, so this function makes the projection to obtain w.

	zrt = octant.depths.get_zrho(Vtransform, Vstretching, km, theta_s, theta_b, 
                    -h_neg_, hc, zeta=zeta, Hscale=3)

	omega[:-1] = omega2w(u_s, v_s, omega, zrt, pm, pn)[:]

	# Change units of w from m/s to m/day

	######### This goes inside omega_to_w
	omega[:] = 24*3600.*omega[:]
	# Cut the values of w up to 15 m/day in absolute value
	omega[np.where(omega > bound)] = bound
	omega[np.where(omega < -bound)] = -bound
	

	# Now interpolate w from sigma to the same z levels of levels_field
	#print '  - Interpolating w to z levels. '

	if z_fields: #This function is not going to be called with z_fields activated, yet.

		w_z = np.zeros((2, ) + s_shape)

		for k in np.arange(roms_levels):
			z = levels_field[k]

			for zrt_k in np.arange(zrt.shape[0] - 1):
				zrtm1 = zrt[zrt_k]
				zrtp1 = zrt[zrt_k + 1]

				in_dz = np.where((zrtm1 < z) & (zrtp1 >= z))

				zrtm1 = zrtm1[in_dz]
				zrtp1 = zrtp1[in_dz]

				if zrtm1.size > 0:
					weigh = abs((z - zrtm1) / (zrtp1 - zrtm1))
					w_z[k][in_dz] = omega[zrt_k][in_dz] + (omega[zrt_k + 1][in_dz] - omega[zrt_k][in_dz]) * weigh
		return w_z
	else:
		return omega


def sigma2z():
	""" 
	Interpolation from sigma to the same z levels of levels_field
	"""
	print '  - Interpolating from sigma to z levels. '


	for i, field in enumerate(fields_s):

		for k, z in enumerate(levels_field):

			for zrt_k in np.arange(roms_levels - 1):

				zrtm1 = zrt[zrt_k]
				zrtp1 = zrt[zrt_k + 1]

				in_dz = np.where((zrtm1 < z) & (zrtp1 >= z))

				zrtm1 = zrtm1[in_dz]
				zrtp1 = zrtp1[in_dz]

				if zrtm1.size > 0:
					weigh = abs((z - zrtm1) / (zrtp1 - zrtm1))
					fields_s[i][k][in_dz] = field[zrt_k][in_dz] + (field[zrt_k + 1][in_dz] - field[zrt_k][in_dz]) * weigh
		#return int_var


def level_limits(levels):
	# Lower and upper values of levels_field if some limit of levels do not belong to levels_field
	levels_wellbounded = levels
	if  not levels[0] in levels_field:
		for j, lev in enumerate(levels_field[:-1]):
			if levels[0] < lev and levels[0] >= levels_field[j + 1]:
				levels_wellbounded.insert(0, levels_field[j])
	if  not levels[-1] in levels_field:
		for j, lev in enumerate(levels_field[:-1]):
			if levels[-1] < lev and levels[-1] >= levels_field[j + 1]:
				levels_wellbounded.append(levels_field[j + 1]) 
	return levels_wellbounded


def cube_marginals(cube, normalize=False, insitu=True):
    """ Taken from http://stackoverflow.com/questions/25494668/best-way-to-plot-a-3d-matrix-in-python and modified"""
    c_fcn = np.mean if normalize else np.sum
    xy = c_fcn(cube, axis=0)
    xz = c_fcn(cube, axis=1)
    yz = c_fcn(cube, axis=2)

    if insitu:
    	xy = cube[0, :, :]
    	xz = cube[:, -1, :]
    	yz = cube[:, :, -1]
    return(xy,xz,yz)


def plotcube(cube_field,x=None,y=None,z=None,normalize=False, insitu=True, plot_front=False,
	 offsets=None, which=(True, True, True), fill=False, bar=False, alpha=1.):

    """
    Use contourf to plot cube marginals
    Taken from http://stackoverflow.com/questions/25494668/best-way-to-plot-a-3d-matrix-in-python and modified
    """

    (Z,Y,X) = cube_field.shape
    (xy,xz,yz) = cube_marginals(cube_field,normalize=normalize, insitu=insitu)
    if x == None: x = np.arange(X)
    if y == None: y = np.arange(Y)
    if z == None: z = np.arange(Z)


    # draw edge marginal surfaces
    if offsets is None:
        offsets = (Z-1,0,X-1) if plot_front else (0, Y-1, 0)

    #print 'cube_marginals shape:', xy.shape ,xz.shape ,yz.shape 
    #print 'offsets:', offsets


    plot_func = ax.contourf if fill else ax.contour

    '''
    cset = plot_func(x[None,:].repeat(Y,axis=0), y[:,None].repeat(X,axis=1), xy, zdir='z', offset=offsets[0], 
                        levels=np.linspace(cube_field.min(), cube_field.max(), 100), cmap=plt.cm.coolwarm, alpha=0.75)
    cset = plot_func(x[None,:].repeat(Z,axis=0), xz, z[:,None].repeat(X,axis=1), zdir='y', offset=offsets[1], 
                        levels=np.linspace(cube_field.min(), cube_field.max(), 100), cmap=plt.cm.coolwarm, alpha=0.75)
    cset = plot_func(yz, y[None,:].repeat(Z,axis=0), z[:,None].repeat(Y,axis=1), zdir='x', offset=offsets[2], 
                        levels=np.linspace(cube_field.min(), cube_field.max(), 100), cmap=plt.cm.coolwarm, alpha=0.75)
    '''


    if fill:
    	vmin = min2
    	vmax = max2
    	#vmin = cube_field.min() if you don't calculate the minima and maxima at the beginning of the script
    	# The same with max and/or mean
    	# The same for the non-filled contours
        if which[2]: cset = ax.contourf(x[None,:].repeat(Y,axis=0), y[:,None].repeat(X,axis=1), xy, zdir='z', offset=offsets[0], 
                            levels=np.linspace(vmin, vmax, nfcontours), cmap=cmap, alpha=alpha,
                             antialiased=True)
        if which[1]: cset = ax.contourf(x[None,:].repeat(Z,axis=0), xz, z[:,None].repeat(X,axis=1), zdir='y', offset=offsets[1], 
                            levels=np.linspace(vmin, vmax, nfcontours), cmap=cmap, alpha=alpha,
                             antialiased=True)
        if which[0]: cset = ax.contourf(yz, y[None,:].repeat(Z,axis=0), z[:,None].repeat(Y,axis=1), zdir='x', offset=offsets[2], 
                            levels=np.linspace(vmin, vmax, nfcontours), cmap=cmap, alpha=alpha,
                             antialiased=True)

        label = labels[fill_field]
        if firstbar_fie and bar:
        	cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        	cbar = fig.colorbar(cset, cax=cbaxes, orientation=cb_orientation)#, fraction = cb_fraction, pad=0.05, panchor=(0., 0.), shrink=cb_shrink)
        	cbar.set_label(label)
        	#ticks = np.linspace(f2plot2.min(), f2plot2.max(), 11) # f2plot2 is the color filled field.
        	ticks = np.linspace(min2, max2, 11)
        	cbar.ax.set_yticks(ticks)
        	cbar.ax.set_yticklabels(['%0.1f'%tick for tick in ticks])
    else:
    	vmin = min1
    	vmax = max1
        if which[2]: cset = ax.contour(x[None,:].repeat(Y,axis=0), y[:,None].repeat(X,axis=1), xy, zdir='z', offset=offsets[0], 
                            levels=np.linspace(vmin, vmax, ncontours), colors=cont_colors, alpha=alpha,
                             linewidth=cont_lw)
        if which[1]: cset = ax.contour(x[None,:].repeat(Z,axis=0), xz, z[:,None].repeat(X,axis=1), zdir='y', offset=offsets[1], 
                            levels=np.linspace(vmin, vmax, ncontours), colors=cont_colors, alpha=alpha,
                             linewidth=cont_lw)
        if which[0]: cset = ax.contour(yz, y[None,:].repeat(Z,axis=0), z[:,None].repeat(Y,axis=1), zdir='x', offset=offsets[2], 
                            levels=np.linspace(vmin, vmax, ncontours), colors=cont_colors, alpha=alpha,
                             linewidth=cont_lw)
    
    # draw wire cube to aid visualization

    if cube:
	    ax.plot([x.max(),x.max(),x.min(),x.min(),x.max()],[y.min(),y.max(),y.max(),y.min(),y.min()],[z.max(),z.max(),z.max(),z.max(),z.max()],
	     c=cube_edg_clrs, lw=4., alpha=cube_edg_alph)
	    ax.plot([x.max(),x.max(),x.max(),x.min(),x.min()],[y.min(),y.min(),y.max(),y.max(),y.max()],[z.max(),z.min(),z.min(),z.min(),z.max()],
	     c=cube_edg_clrs, lw=4., alpha=cube_edg_alph)
	    ax.plot([x.max(),x.max()],[y.max(),y.max()],[z.min(),z.max()], c=cube_edg_clrs, lw=4.)

	    npoints = 800

	    aristaz_x = np.ones(npoints) * x.max()
	    aristaz_y = np.ones(npoints) * y.max()
	    aristaz_z = np.linspace(z.min(), z.max(), npoints)
	    ax.scatter(aristaz_x, aristaz_y, aristaz_z, marker='o', s=3., color=cube_edg_clrs, linewidth=0., alpha=cube_edg_alph)

	    aristax_x = np.linspace(x.min(), x.max(), npoints)
	    aristax_y = np.ones(npoints) * y.max()
	    aristax_z = np.ones(npoints) * z.max()
	    ax.scatter(aristax_x, aristax_y, aristax_z, marker='o', s=3., color=cube_edg_clrs, linewidth=0., alpha=cube_edg_alph)

	    aristay_x = np.ones(npoints) * x.max()
	    aristay_y = np.linspace(y.min(), y.max(), npoints)
	    aristay_z = np.ones(npoints) * z.max()
	    ax.scatter(aristay_x, aristay_y, aristay_z, marker='o', s=3., color=cube_edg_clrs, linewidth=0., alpha=cube_edg_alph)


def new_grid(son, sox, san, sax, x, y, field_list, new_shape=(100, 100), only_xy=False):

	"""
	"""

	x_rect, y_rect = np.meshgrid(np.linspace(son, sox, new_shape[1]), \
	                            np.linspace(san, sax, new_shape[0]))
	#print x_rect.shape
	#print x_rect[0]
	#print x_rect[:, 0]
	if only_xy:
		return x_rect, y_rect
	else:
		coord_full = np.array([x.flatten(),y.flatten()])
		coord0 = np.transpose(np.array([x_rect.flatten(),y_rect.flatten()]))
		fields_rect = []

		for field in field_list:
			print field.shape

			if field.shape == z_shape_u:
				coord = np.transpose(np.array([x[:, :-1].flatten(),y[:, :-1].flatten()]))
			elif field.shape == z_shape_v:
				coord = np.transpose(np.array([x[:-1, :].flatten(),y[:-1, :].flatten()]))
			else:
				coord = np.transpose(coord_full)

			if field.ndim == 2:
				tmp = np.zeros(new_shape)
				tmp[:] = griddata(coord,field.flatten(),coord0).reshape(new_shape)
			elif field.ndim == 3: # Be careful! axis=0 could be either time (for zeta) or sigma level (for zrt), depending on the variable.
				tmp = np.zeros(field.shape[:1] + new_shape)
				#print coord.shape, field[0].flatten().shape, field[0].shape, coord0.shape
				for t in np.arange(field.shape[0]):
					tmp[t] = griddata(coord,field[t].flatten(),coord0).reshape(new_shape)
				#fields_rect.append(tmp)
			elif field.ndim == 4:
				tmp = np.zeros(field.shape[:2] + new_shape)
				for t in np.arange(field.shape[0]):
					for lev in np.arange(field.shape[1]):
						#print t, lev, coord.shape, field[t, lev].flatten().shape, coord0.shape
						tmp[t, lev] = griddata(coord,field[t, lev].flatten(),coord0).reshape(new_shape)
			## Check for nan's, mask them and turn them into zero.
			tmp = np.ma.array(tmp, mask=np.ma.masked_where(tmp, np.isnan(tmp)))
			tmp.data[np.where(np.isnan(tmp.data))] = 0.
			fields_rect.append(tmp)

		return x_rect, y_rect, fields_rect


def get_region(son, sox, san, sax):

	# Corners of the domain to be plotted
	nw = (son, sax)
	ne = (sox, sax)
	se = (sox, san)
	sw = (son, san)

	corners = [nw, ne, se, sw]

	icorner = np.zeros(4, dtype=int)
	jcorner = np.zeros(4, dtype=int)

	i_surr4 = np.zeros((4,4), dtype=int)
	j_surr4 = np.zeros((4,4), dtype=int)

	# Find the closest point in the grid to each corner
	for k, corner in enumerate(corners):

		icorner[k], jcorner[k], _ = nearest_2D_single(corner[0], corner[1], _lon, _lat)
		i_surr4[k, :], j_surr4[k, :] = surrounding4(corner[0], corner[1], _lon, _lat)

	# Eliminate discrepancies between the chosen points
	'''minlonindex = icorner.min()
	maxlonindex = icorner.max()
	minlatindex = jcorner.min()
	maxlatindex = jcorner.max()'''

	minlonindex = i_surr4.min()
	maxlonindex = i_surr4.max()
	minlatindex = j_surr4.min()
	maxlatindex = j_surr4.max()

	# Finally, set the region
	# This region MUST contain the four corners and in fact, it does
	region = np.s_[minlonindex:maxlonindex + 1, minlatindex:maxlatindex + 1]
	region_u = np.s_[minlonindex:maxlonindex + 1, minlatindex:maxlatindex]
	region_v = np.s_[minlonindex:maxlonindex, minlatindex:maxlatindex + 1]

	return region, region_u, region_v


def nearest(lon_pt, lat_pt, lon2d, lat2d, four=False, test=False):
    """
    Return the nearest i, j point to a given lon, lat point
    in a ROMS grid
        i,j = nearest(lon_pt,lat_pt,lon2d,lat2d,four=False,test=False)
    If four is not None then the four nearest (i, j)
    points surrounding the point of interest are returned
    To test with a figure put test=True
    By Evan Mason
    """
    def test_it(i, j, lon_pt, lat_pt, lon2d, lat2d):
        print '---Testing the locations...'
        i1, i2, i3, i4 = i
        j1, j2, j3, j4 = j
        '''plt.figure(1)
        plt.plot(lon2d, lat2d, 'ro')
        plt.text(lon_pt, lat_pt, '0')
        plt.axis('image')
        plt.text(lon2d[j1, i1], lat2d[j1, i1], 'A')
        plt.text(lon2d[j2, i2], lat2d[j2, i2], 'B')
        plt.text(lon2d[j3, i3], lat2d[j3, i3], 'C')
        plt.text(lon2d[j4, i4], lat2d[j4, i4], 'D')
        plt.title('Is the 0 in the middle?')
        plt.show()'''
        return
    def four_points(i1, j1, lon_pt, lat_pt, lon2d, lat2d):
        iii   = np.array([1, -1, -1,  1])
        jjj   = np.array([1,  1, -1, -1])
        big_dist = 1e99 # arbitrary large distance
        for ii, jj in zip(iii, jjj):
            dist2 = np.hypot(lon2d[j1 + jj, i1 + ii] - lon_pt,
                             lat2d[j1 + jj, i1 + ii] - lat_pt)
            if dist2<big_dist:
                big_dist = dist2
                i2, j2 = i1+ii, j1+jj
        if i2>i1: i3 = i2; i4 = i1
        else:     i3 = i2; i4 = i1
        if j2>j1: j3 = j1; j4 = j2
        else:     j3 = j1; j4 = j2
        return np.array([i1, i2, i3, i4]), \
               np.array([j1, j2, j3, j4])
    # Main...
    d    = np.hypot(lon2d-lon_pt, lat2d-lat_pt)
    j, i = np.unravel_index(d.argmin(), d.shape)
    dist = d[j,i] # j,i not i,j
    
    if four is True:
        #print '---Getting 4 surrounding points...'
        i, j = four_points(i, j, lon_pt, lat_pt, lon2d, lat2d)
        if test is True: test_it(i, j, lon_pt, lat_pt, lon2d, lat2d)
        return i, j
    else:
        return i, j


def create_slices(field1, field2):
	# Create the slices, interpolating if necessary
	ii = 0
	jj = 0
	for i, level in enumerate(levels):			
		if level in levels_field:
			#print 'level in levels_field', i, ii, jj, level  				
			field_slice[i] = field1[ii]
			field_slice2[i] = field2[ii]
			ii = ii + 1
			jj = 0
		elif ii != 0:
			#print 'level not in levels_field', i, ii, jj, level
			jj = jj + 1
			weight = (level - levels[old_levels_inds[ii - 1]])/(levels[old_levels_inds[ii]] - levels[old_levels_inds[ii - 1]])
			field_slice[i, :] = (field1[ii - 1]) + (field1[ii] - field1[ii - 1])*weight
			field_slice2[i, :] = (field2[ii - 1]) + (field2[ii] - field2[ii - 1])*weight
		else:
			print 'ERROR: First level outside levels_field.', ii, level
			break

	#return field_slice, field_slice2


def wall_slices(which):
	if moving_v_slice:
		if which_slices[0]: offsets = (lfc[0], y_eddy_r.max(), wall_offset)
		if which_slices[1]: offsets = (lfc[0], wall_offset, x_eddy_r.max())
	else:
		offsets = (lfc[0], y_eddy_r.max(), x_eddy_r.max())

	plotcube(f2plot1_int,x=x_eddy_r[0,:],y=y_eddy_r[:,0],z=lfc,normalize=True, insitu=True, plot_front=True, 
			offsets=offsets, which=which, fill=False, alpha=alpha_slices)
	plotcube(f2plot2_int,x=x_eddy_r[0,:],y=y_eddy_r[:,0],z=lfc,normalize=True, insitu=True, plot_front=True, 
			offsets=offsets, which=which, fill=True, bar=show_bar, alpha=alpha_slices)


def slice_movement(t, mov_type):
	if 'Sin' in mov_type:
		return np.abs(np.sin((t/ntimes) * 1. * np.pi))
	if 'Lin' in mov_type:
		return 1. - 2. * np.abs((t / ntimes - 0.5))
	if 'Both' in mov_type:
		if t < t4lon_slice:
			return 1. - 2. * np.abs((2. * t / ntimes - 0.25))
		else:
			return 1. - 2. * np.abs(((2. * t - t4lon_slice) / ntimes - 0.25))


def draw_fields(f1, f2, xc=None, yc=None, v_raw=None, what2draw=None, offsets=None):

	if what2draw == 'Slices':
		if show_black_field:
			# v_raw[0] if you use show_slice_imit_isosurfs instead of show_slice_imit_isosurfs2
			show_slice_imit_isosurfs2(x_eddy_r, y_eddy_r, f1, min1, max1, level, firstbar_fie, cmap=None, 
				alpha=alpha_slices, label=labels[cont_field], fill=False, bar=show_bar_cont, colors=cont_colors, zorder=2)
		if show_colored_field:
			# v_raw[1] if you use show_slice_imit_isosurfs instead of show_slice_imit_isosurfs2
			show_slice_imit_isosurfs2(x_eddy_r, y_eddy_r, f2, min2, max2, level, firstbar_fie, cmap=cmap, 
				alpha=alpha_slices, label=labels[fill_field], fill=fill_contours, bar=show_bar, pos_neg=False, zorder=2)
	elif what2draw == 'Cube':
		if show_black_field:
			plotcube(f1,x=xc,y=yc,z=lfc,normalize=True, insitu=True, plot_front=True, 
					offsets=offsets, fill=False, alpha=alpha_slices)
		if show_colored_field:
			plotcube(f2,x=xc,y=yc,z=lfc,normalize=True, insitu=True, plot_front=True, 
					offsets=offsets, fill=fill_contours, bar=show_bar, alpha=alpha_slices)


#####################################################################
### Control Pannel. ###
#####################################################################

##	Select the region we want to plot

bathymetry = True
bath_res_fac = 30

# Limit values for the basin
son_basin = -20. #87.	#-20.	#-16.738680588838935
sox_basin = -12. #99.	#-13.	#-15.963741383221457
san_basin = 24.02 #22.	#25. 	#26.049154597891373
sax_basin = 31. #31.	#30. 	#26.589057212058627

# Limit values for the tracks (for the eddy "jail"). These values are only going to be used if not custom_region

if bathymetry:
	'''son_eddy = son_basin
	sox_eddy = sox_basin
	san_eddy = san_basin
	sax_eddy = sax_basin'''
	son_eddy = -16.838680588838935 #87.	#-20.	#-16.738680588838935
	sox_eddy = -15.563741383221457 #99.	#-13.	#-15.963741383221457
	san_eddy = 25.949154597891373 #22.	#25. 	#26.049154597891373
	sax_eddy = 26.989057212058627 #31.	#30. 	#26.589057212058627
else:
	son_eddy = -16.838680588838935 #87.	#-20.	#-16.738680588838935
	sox_eddy = -15.563741383221457 #99.	#-13.	#-15.963741383221457
	san_eddy = 25.949154597891373 #22.	#25. 	#26.049154597891373
	sax_eddy = 26.989057212058627 #31.	#30. 	#26.589057212058627

### Other parameters

matplotlib.rcParams.update({'font.size':18,'font.family':'STIXGeneral','mathtext.fontset':'stix','text.usetex':False})
dpi = 100
figsize = (19.2, 10.8)

fieldtimes = 2 # How many ROMS outputs. Please, two or more.
z_fields = False
fls_from_tracpy = True
less_particles = False
less_factor = 2
reduce_w = False
custom_region = True

show_axes = True
show_time = False

show_particles = True
parts_color, parts_edge_col, patrs_edge_size = [0.4, 0.4, 0.4], [0.2, 0.2, 0.2], 0.2
parts_size = 2.
initial_subregion = False
fade_initial = False
fade_outs_custom_reg = True
fade_field = None # Make it work for different fields
want_subregion = False
# subregion is set in line 682.

# Fields that we want to use
want = {'w':False, 'uv':False, 'h_vorticity':False, 'pt':True, 'salt':False, 'Density':False}
labels = {'w':'w (m/s)', 'uv':'u (m/s)', 'h_vorticity':'Vorticity (s-1)', 'pt':'Temperature (Celsius)', 
'salt':'Salinity (g/kg)', 'Density':'Density (kg/m3)'}
cont_field = 'pt'
fill_field = 'pt'

fields_calc = True
fields_plot = True
show_bar = True
show_bar_cont = False
cb_orientation = "vertical"
cb_fraction = 0.03
cb_shrink = 0.7
cb_ticks = 5
cmap = cm.jet
#cmap.set_bad('magenta')
cont_colors = 'k' # None: if you want a cmap, or a certain color for that color.
cont_lw = 0.5
ncontours = 30
nfcontours = 200
show_black_field = False
show_colored_field = True

fade_method = None #'Masking'
colorby = None #'fixed' #'field'
fvalue = 1026.98
fwidth = 0.04

zoom = False
move_cam = True
height = 100
azimut = 360

## Levels for the slices. Too complicated this way...
levels = [-100, -130, -150] 	# Levels of the slices that I want to plot. The first 
									# and the last must be in levels_field and are not plotted.
z_levels = 2 # Number of z levels that we want to obtain from the sigma levels.
cube = False
cube_edg_clrs, cube_edg_alph = 'k', 1.
horiz_slices = True
vertic_slices = False
which_slices = (True, False, False) # (zdir:longitude, zdir:latitude, zdir:vertical)
moving_h_slice = False
moving_v_slice = False
mov_type = 'Linear'
fill_contours = False

delta_lon = sox_eddy - son_eddy 	# No need to change these two lines; they are just calculations.
delta_lat = sax_eddy - san_eddy

mhs_origin = -150
mhs_finish = -100

mvs_origin_lon = sox_eddy #- delta_lon * 0.5
mvs_finish_lon = son_eddy #- delta_lon * 0.5
if which_slices[0] and not which_slices[1]: wall_offset = mvs_origin_lon

mvs_origin_lat = sax_eddy
mvs_finish_lat = san_eddy
if which_slices[1] and not which_slices[0]: wall_offset = mvs_origin_lat

fade_slices = False
alpha_slices = 1.

# Limits of the figure
z_lims = (10, -5000)
x_lims = (son_basin, sox_basin)
y_lims = (san_basin, sax_basin)

framesperday = 48
anim_nu = 157

t0 = 120
tsep = 30


#####################################################################
### Create a folder for the movie ###
#####################################################################

anim_name = 'Animation_%04i'%anim_nu
folder_path = '/home/mcapllonch/Documents/Animacions/' + anim_name + '/'

if not os.path.exists(folder_path):
	os.makedirs(folder_path)

#####################################################################
### Names of the files ###
#####################################################################

dslice = slice(0, 60, 1)
fslice = slice(0, 10, 10)

## Names of the files that have the data:
	# Floats
fl_evan = '/home/mcapllonch/Downloads/python_scripts/uswc_floats_PUMP8_center.nc'
fl_evan2 = '/predator/emason/roff/na75/avgs_3day/averages2_bwds_feb.nc'
fl_tracpy = '/home/mcapllonch/scripts_miquel/tracks/prova_24.nc'
	# Fields
fi_z_roms = '/predator/emason/runs2009/gc_2009_1km_60/z_roms_avg.0990.nc'
fi_roms = '/predator/emason/runs2009/gc_2009_1km_60/roms_avg.0990.nc'
fi_roms2 = '/predator/emason/runs2009/na_2009_7pt5km/roms_avg.4080.nc'
fi_roms3 = '/home/mcapllonch/scripts_miquel/nc_files_for_tracpy/roms_fileds.nc'
	# Grids
grd_roms = '/predator/emason/runs2009/gc_2009_1km_60/gc_2009_1km_grd_smooth.nc'
grd_roms2 = '/predator/emason/runs2009/na_2009_7pt5km/roms_grd_NA2009_7pt5km.nc'
grd_roms3 = '/home/mcapllonch/scripts_miquel/nc_files_for_tracpy/roms_grid.nc'



#################################################################
### Grid, shape of the variables, ocean time and vertical z levels ###
#################################################################

print '*** Reading grid, shape of the variables, ocean time and vertical z levels. ***'
print

# Read grid variables
with Dataset(grd_roms3) as nc:
	_lon = nc.variables['lon_rho'][:]
	_lat = nc.variables['lat_rho'][:]
	h_neg = -nc.variables['h'][:]
	f = nc.variables['f'][:]
	pm = nc.variables['pm'][:]
	pn = nc.variables['pn'][:]
	mask_rho = nc.variables['mask_rho'][:]

# Get the shape of a variable from ROMS and some relevant parameters
roms_file = fi_z_roms if z_fields else fi_roms3

with Dataset(roms_file) as nc:
	roms_times = nc.variables['ocean_time'][:].size
	if not fls_from_tracpy: ocean_time = nc.variables['ocean_time'][:]
	if z_fields:
		roms_levels = nc.dimensions['depth']
	else:
		if 's_w' in nc.variables:
		    sc_r = nc.variables['s_w'][:]
		else:
		    sc_r = nc.variables['sc_w'][:]
		roms_levels = sc_r.size
	theta_s = nc.variables['theta_s'][:]
	theta_b = nc.variables['theta_b'][:]
	hc = nc.variables['hc'][:]

# Some parameters

km = sc_r.shape[0] - 1 
Vtransform = 2
Vstretching = 2

# Levels
if z_fields:
	# Levels of the z_roms output
	levels_field = [0, -5, -10, -15, -20, -30, -40, -50, -75, -100, -200, -400, 
		-600, -800, -1000, -1250, -1500, -1750, -2000, -2250, -2500, -2750, -3000]
else:
	levels_field = np.linspace(levels[0], levels[-1], z_levels).tolist()

# Lower and upper values of levels_field if some limit of levels do not belong to levels_field

#levels = level_limits(levels)

if  not levels[0] in levels_field:
	for j, lev in enumerate(levels_field[:-1]):
		if levels[0] < lev and levels[0] >= levels_field[j + 1]:
			levels.insert(0, levels_field[j])
if  not levels[-1] in levels_field:
	for j, lev in enumerate(levels_field[:-1]):
		if levels[-1] < lev and levels[-1] >= levels_field[j + 1]:
			levels.append(levels_field[j + 1]) 

#################################################################
### Particles ###
#################################################################

print '*** Reading particles. ***'
print

# Read
if fls_from_tracpy:
	with Dataset(fl_tracpy) as nc:
		lats_nc = nc.variables['latp'][:].T
		lons_nc = nc.variables['lonp'][:].T
		depths_nc = nc.variables['zp'][:].T
else:
	with Dataset(fl_evan) as nc:
		lats_nc = nc.variables['lat'][:]
		lons_nc = nc.variables['lon'][:]
		depths_nc = nc.variables['depth'][:]	
		ocean_time_floats = nc.variables['ocean_time'][dslice]

### Select only some of the particles if there are too many

if less_particles:
	lats_nc = lats_nc[:, ::less_factor]
	lons_nc = lons_nc[:, ::less_factor]
	depths_nc = depths_nc[:, ::less_factor]


### Removing bad values from the particles

print '  - Removing bad values of the particles. '

## Remove nans and turn them into 10000. to make the arrays look 
## like they were nottaken from tracpy

inds_o = np.where(np.isnan(lons_nc))
inds_d = np.where(np.isnan(depths_nc))

lons_nc[inds_o] = 10000.
lats_nc[inds_o] = 10000.
depths_nc[inds_o] = 10000.
lons_nc[inds_d] = 10000.
lats_nc[inds_d] = 10000.
depths_nc[inds_d] = 10000.

## Creating the mask for the bad values of the arrays:
bad_lons = np.ma.masked_greater(lons_nc, 1000.).mask
bad_lats = np.ma.masked_greater(lats_nc, 1000.).mask
bad_depths = np.ma.masked_greater(depths_nc, 1000.).mask
badvalues_mask0 = np.ma.mask_or(bad_lons, bad_lats)
badvalues_mask = np.ma.mask_or(badvalues_mask0, bad_depths)

## Check and select the coincidence of ocean times between the floats and the
## fields
if fls_from_tracpy:
	combined_mask = badvalues_mask
else:
	i, j = np.meshgrid(*map(np.arange, lats_nc.shape), indexing='ij')
	days_mask = (i >= roms_times) 	# That must be checked with ocean_times
	combined_mask = np.ma.mask_or(days_mask, badvalues_mask)

lats_nc = np.ma.array(lats_nc, mask=combined_mask)
lons_nc = np.ma.array(lons_nc, mask=combined_mask)
depths_nc = np.ma.array(depths_nc, mask=combined_mask)

# Reduce vertical displacements if they are too high. Not necesary anymore if tracpy works fine.
if reduce_w:

	print '*** Reduction of the vertical velocities of the particles. ***'
	print

	dd = depths_nc[1:, :] - depths_nc[:-1, :]
	maxdd = np.array([abs(dd.min()), abs(dd.max())]).max()
	depths_nc = no_excursion(framesperday, maxdd)


### Masking particles outside a certain subregion ###

## Choose ALL the floats (that's why to use _nc) that are inside the region:
## Unmask the arrays
print '  - Masking the floats outside a certain subregion.'
# These arrays must be unmasked. Check if you can avoid masking them before.
lons_nc = lons_nc.data
lats_nc = lats_nc.data
depths_nc = depths_nc.data

## Create a new combined mask and mask the unmasked arrays
## to restrict the particles to those inside the domain we want to plot

subregion = np.array([[son_eddy, sox_eddy], [san_eddy, san_eddy], [-125., -120.]])
try:
	if want_subregion:
		lons_mask = np.ma.masked_outside(lons_nc, subregion[0, 0], subregion[0, 1]).mask
		lats_mask = np.ma.masked_outside(lats_nc, subregion[1, 0], subregion[1, 1]).mask
		depths_mask = np.ma.masked_outside(depths_nc, subregion[2, 0], subregion[2, 1]).mask
	else:
		lons_mask = np.ma.masked_outside(lons_nc, son_eddy, sox_eddy).mask
		lats_mask = np.ma.masked_outside(lats_nc, san_eddy, sax_eddy).mask
		depths_mask = np.ma.masked_outside(depths_nc, z_lims[-1], z_lims[0]).mask
except ValueError:
	print('ERROR: "Corners do not fit in the model domain" is no longer what is happenning here.')
	sys.exit()
combined_mask_2 = np.ma.mask_or(lons_mask, lats_mask)
combined_mask_a = np.ma.mask_or(depths_mask, depths_mask)
combined_mask_3 = np.ma.mask_or(combined_mask_a, combined_mask_2)
combined_mask_4 = np.ma.mask_or(combined_mask, combined_mask_3)
lons_domain = np.ma.array(lons_nc, mask=combined_mask_4, hard_mask=True)
lats_domain = np.ma.array(lats_nc, mask=combined_mask_4, hard_mask=True)
depths_domain = np.ma.array(depths_nc, mask=combined_mask_4, hard_mask=True)

# If we want a subregion only at the beginning:

if initial_subregion:
	lons_mask = np.ma.masked_outside(lons_nc[1, :], subregion[0, 0], subregion[0, 1]).mask
	lats_mask = np.ma.masked_outside(lats_nc[1, :], subregion[1, 0], subregion[1, 1]).mask
	depths_mask = np.ma.masked_outside(depths_nc[1, :], subregion[2, 0], subregion[2, 1]).mask
	combined_mask_2 = np.ma.mask_or(lons_mask, lats_mask)
	combined_mask_a = np.ma.mask_or(depths_mask, depths_mask)
	combined_mask_3 = np.ma.mask_or(combined_mask_a, combined_mask_2)
	combined_mask_4 = np.ma.mask_or(combined_mask[1, :], combined_mask_3)

	lons_good = np.ma.array(lons_nc[1, :], mask=combined_mask_4, hard_mask=True)

	lg = lons_good.compressed()
	inds = np.asarray([np.where(lons_domain[1, :] == lg[i])[0] for i in range(lg.size)]).flatten()
	lons_domain = lons_domain[:, inds]
	lats_domain = lats_domain[:, inds]
	depths_domain = depths_domain[:, inds]

print


###############################################
###	Set the region we want to plot ###
###############################################

print '*** Setting the region we want to plot. ***'
print

if not custom_region:
	son_eddy = lons_nc.min()
	sox_eddy = lons_nc.max()
	san_eddy = lats_nc.min()
	sax_eddy = lats_nc.max()		

region_eddy, rg_e_u, rg_e_v = get_region(son_eddy, sox_eddy, san_eddy, sax_eddy)
x_eddy, y_eddy = _lon[region_eddy], _lat[region_eddy]

dmax = np.sqrt(np.asarray(np.gradient(y_eddy)).max()**2 + np.asarray(np.gradient(x_eddy)).max()**2)

if bathymetry:
	region_basin, rg_b_u, rg_b_v = get_region(son_basin, sox_basin, san_basin, sax_basin)
	x_basin, y_basin = _lon[region_basin], _lat[region_basin]
	h_neg_basin = h_neg[region_basin]

# Indices for the vertical levels to set the region
lfil0 = levels_field.index(levels[0])
lastlev = levels[-1]
for i, lev in enumerate(levels_field[:-1]):
	if lastlev < lev and lastlev >= levels_field[i + 1]:
		lastlev = levels_field[i + 1]

lfilm1 = levels_field.index(lastlev)
lev_fie_chosen = levels_field[lfil0:lfilm1+1]
lfc = np.asarray(lev_fie_chosen)

# Region that includes a depth range (between the slices, now)
# and the horizontal region created before

vert_slice = np.s_[lfil0:lfilm1+1] if z_fields else np.s_[:]
region3d_eddy = (vert_slice, ) + region_eddy 
rg3d_u_eddy = (vert_slice, ) + rg_e_u
rg3d_v_eddy = (vert_slice, ) + rg_e_v


###############################################
###	Define variables and allocate memory ###
###############################################

print '*** Define variables and allocate memory. ***'
print

# First of all, cut some grid variables to the region
# I don't read them like this from the grid file for in case I want to plot/use them in a wider area.
f = f[region_eddy]
pm = pm[region_eddy]
pn = pn[region_eddy]
mask_rho = mask_rho[region_eddy]

# Now we can proceed...

s_shape = (2, km, ) + x_eddy.shape
s_shape_u = (2, km, ) + (x_eddy.shape[0], x_eddy.shape[1] - 1) 
s_shape_v = (2, km, ) + (x_eddy.shape[0] - 1, x_eddy.shape[1]) 
s_shape_w = (2, roms_levels, ) + x_eddy.shape
z_shape = (2, z_levels, ) + x_eddy.shape
z_shape_u = (2, z_levels, ) + (x_eddy.shape[0], x_eddy.shape[1] - 1)
z_shape_v = (2, z_levels, ) + (x_eddy.shape[0] - 1, x_eddy.shape[1])
if moving_v_slice:
	if which_slices[0] and not which_slices[1]:
		new_shape = (100, 1000)
	elif which_slices[1] and not which_slices[0]:
		new_shape = (1000, 100)
	elif which_slices[0] and which_slices[1]:
		new_shape = (1000, 1000)

else:
	new_shape = (100, 100)
rect_shape = (2, z_levels, ) + new_shape

#x_eddy_r = np.zeros(new_shape)
#y_eddy_r = np.zeros(new_shape)

x_eddy_r, y_eddy_r, lmask = new_grid(son_eddy, sox_eddy, san_eddy, sax_eddy, 
					x_eddy, y_eddy, [mask_rho], new_shape=new_shape, only_xy=False)
thereis_land = np.in1d(True, mask_rho.flatten())[0]
mask_newshape = lmask[0]
mask_newshape[np.where(mask_newshape > 0.)] = 1.

mask_newgrid = np.empty(rect_shape)
for time in np.arange(2):
	for lev in np.arange(z_levels):
		mask_newgrid[time, lev] = mask_newshape[:]

if fields_calc:
	# Fields that are going to be used at some point. The list is for setting up the variables needed for the interpolation to a new grid.

	fields_s = []
	fields_z = []

	zeta = np.array([np.empty_like(x_eddy), np.empty_like(x_eddy)])
	fields_z.append(zeta)
	zeta_rect = np.zeros((2, ) + new_shape)
	zrt = np.zeros(s_shape)

	if want['w']:
		omega = np.zeros(s_shape_w)
		w = np.zeros(s_shape_w)
		fields_s = fields_s + [w]
		omega_z = np.zeros(z_shape)
		w_z = np.zeros(z_shape)
		fields_z = fields_z + [w_z]
		omega_r = np.zeros(rect_shape)
		w_r = np.ma.array(np.zeros(rect_shape), mask=mask_newgrid)
	if want['w'] or want['uv'] or want['h_vorticity']:
		'''u = np.zeros(s_shape_u)
		v = np.zeros(s_shape_v)'''
		u = np.zeros(s_shape)
		v = np.zeros(s_shape)
		'''u_z = np.zeros(z_shape_u)
		v_z = np.zeros(z_shape_v)'''
		u_z = np.zeros(z_shape)
		v_z = np.zeros(z_shape)
		fields_s = fields_s + [u, v]
		fields_z = fields_z + [u_z, v_z]
		u_r = np.ma.array(np.zeros(rect_shape), mask=mask_newgrid)
		v_r = np.ma.array(np.zeros(rect_shape), mask=mask_newgrid)
	if want['h_vorticity']:
		hvort = np.empty_like(u_r)
	if want['pt'] or want['Density']:
		pt = np.zeros(s_shape)
		pt_z = np.zeros(z_shape)
		pt_r = np.ma.array(np.zeros(rect_shape), mask=mask_newgrid)
		fields_s = fields_s + [pt]
		fields_z = fields_z + [pt_z]
		combined_mask_6 = np.empty_like(pt_r)
	if want['salt'] or want['Density']:
		salt = np.zeros(s_shape)
		salt_z = np.zeros(z_shape)
		salt_r = np.ma.array(np.zeros(rect_shape), mask=mask_newgrid)
		fields_s = fields_s + [salt]
		fields_z = fields_z + [salt_z]
		combined_mask_6 = np.empty_like(salt_r)
	if want['Density']:
		press = np.zeros(s_shape[1:]) 	# Only a tmp variable
		press_r = np.zeros(pt_r.shape[1:]) 	# Only a tmp variable
		dens = np.zeros(s_shape)
		dens_r = np.ma.array(np.zeros(rect_shape), mask=mask_newgrid)

	# Set the fields to plot and slices:

	f2plot1 = np.zeros(rect_shape)
	f2plot2 = np.zeros(rect_shape)

	f2plot1_int = np.zeros(rect_shape[1:])
	f2plot2_int = np.zeros(rect_shape[1:])

	field_slice = np.zeros((len(levels), ) + new_shape)
	field_slice2 = np.zeros((len(levels), ) + new_shape)

	## Set the indices for the regions of the z levels that are over topography now if zrt_indep_t = True

	zrt_indep_t = False

	if zrt_indep_t:
		# CAREFUL! This is not going to work.
		zrt[:] = octant.depths.get_zrho(Vtransform, Vstretching, km, theta_s, theta_b, 
		                    -h_neg[region_eddy], hc, zeta=0., Hscale=3)[:]

		zrtm1 = range(len(levels_field))
		zrtp1 = range(len(levels_field))
		weigh = range(len(levels_field))
		in_dz = range(len(levels_field))

		for kz, z in enumerate(levels_field):
			for ks in np.arange(zrt.shape[0] - 1):

				in_dz[kz] = np.where((zrt[ks] < z) & (zrt[ks + 1] >= z))

				zrtm1[kz] = zrt[ks][in_dz[kz]]
				zrtp1[kz] = zrt[ks + 1][in_dz[kz]]

				if zrtm1[kz].size > 0:
					weigh[kz] = abs((z - zrtm1[kz]) / (zrtp1[kz] - zrtm1[kz]))


###############################################
###	Means, max and mis ###
###############################################

if fields_calc:
	print 'Calculating means, maxima and minima'

	means = []
	maxims = []
	minims = []

	if want['w']:
		w_means = np.zeros(roms_times)
		w_maxims = np.zeros(roms_times)
		w_minims = np.zeros(roms_times)
	if want['w'] or want['uv'] or want['h_vorticity']:
		u_means = np.zeros(roms_times)
		u_maxims = np.zeros(roms_times)
		u_minims = np.zeros(roms_times)
		v_means = np.zeros(roms_times)
		v_maxims = np.zeros(roms_times)
		v_minims = np.zeros(roms_times)
	if want['pt'] or want['Density']:
		pt_means = np.zeros(roms_times)
		pt_maxims = np.zeros(roms_times)
		pt_minims = np.zeros(roms_times)
	if want['salt'] or want['Density']:
		salt_means = np.zeros(roms_times)
		salt_maxims = np.zeros(roms_times)
		salt_minims = np.zeros(roms_times)

	for t_r in np.arange(roms_times):
		print t_r, ' of ', roms_times	
		region4d_eddy = (t_r, ) + region3d_eddy
		## Read Fields
		with Dataset(roms_file) as nc:
			zeta[0] = nc.variables['zeta'][t_r][region_eddy]
			zrt[0] = octant.depths.get_zrho(Vtransform, Vstretching, km, theta_s, theta_b, 
				                    -h_neg[region_eddy], hc, zeta=zeta[0], Hscale=3)[:]
			inds = np.where((zrt[0] < levels[0]) & (zrt[0] > levels[-1]))
			if want['pt'] or want['Density']:
				pt[0] = nc.variables['temp'][region4d_eddy]
				pt_means[t_r] = pt[0][inds].mean()
				pt_maxims[t_r] = pt[0][inds].max()
				pt_minims[t_r] = pt[0][inds].min()
			if want['salt'] or want['Density']:	
				salt[0] = nc.variables['salt'][region4d_eddy]
				salt_means[t_r] = salt[0][inds].mean()
				salt_maxims[t_r] = salt[0][inds].max()
				salt_minims[t_r] = salt[0][inds].min()
			if want['w'] or want['uv'] or want['h_vorticity']:
				u[0, :, :, :-1] = nc.variables['u'][region4d_eddy]
				u_means[t_r] = u[0][inds].mean()
				u_maxims[t_r] = u[0][inds].max()
				u_minims[t_r] = u[0][inds].min()
				v[0, :, :-1, :] = nc.variables['v'][region4d_eddy]
				v_means[t_r] = v[0][inds].mean()
				v_maxims[t_r] = v[0][inds].max()
				v_minims[t_r] = v[0][inds].min()	

	if want['w'] or want['uv'] or want['h_vorticity']:
		u_mean = u_means.mean()
		u_maxim = u_maxims.max()
		u_minim = u_minims.min()
		v_mean = v_means.mean()
		v_maxim = v_maxims.max()
		v_minim = v_minims.min()	
	if want['pt'] or want['Density']:
		pt_mean = pt_means.mean()
		pt_maxim = pt_maxims.max()
		pt_minim = pt_minims.min()
	if want['salt'] or want['Density']:	
		salt_mean = salt_means.mean()
		salt_maxim = salt_maxims.max()
		salt_minim = salt_minims.min()


	## Auxiliary variable:
	abv = np.zeros(2)

#####################################################################################
### Choose if we want to mask particles outside a field range at the initial time ###
#####################################################################################

#print 'Choosing if we want to mask particles outside a field range at the initial time:', fade_initial

if fade_initial:

	with Dataset(roms_file) as nc:
		if want['pt'] or want['Density']:
			pt[0] = nc.variables['temp'][t_roms][region3d_eddy]
		if want['salt'] or want['Density']:
			salt[0] = nc.variables['salt'][t_roms][region3d_eddy]
		## You should calculate density here.

	_, lons_good = fade(x_eddy, y_eddy, lons_domain[1].compressed(), lats_domain[1].compressed(), 
		depths_domain[1].compressed(), 	fvalue = fvalue, fwidth = fwidth, 
		fade_field=dens[0], fade_method='Masking')
	lg = lons_good.compressed()
	inds = np.asarray([np.where(lons_domain[1, :] == lg[i])[0] for i in range(lg.size)]).flatten()
	if inds.dtype in 'float64':
		inds.dtype = 'int64'
	lons_domain = lons_domain[:, inds]
	lats_domain = lats_domain[:, inds]
	depths_domain = depths_domain[:, inds]

###############################################
###	MASTER TIME LOOP ###
###############################################

print '### MASTER TIME LOOP ###'

## Some stepping parameters
if fls_from_tracpy:
	ntimes = np.array([lons_domain.shape[0], (roms_times - 1) * framesperday]).min()
	print 'Floats from tracpy. Number of time samples for them:', ntimes
else:
	ntimes = ocean_time.size

nsteps = ntimes // (roms_times - 1) # Number of particle steps per ROMS step. 
nisteps = nsteps # Number of interpolated ROMS steps per ROMS step. Not necessarily equal to nsteps, because this refers to the frames.
nnew = nisteps

# Some parameters for the animation

# Camera
height_cam = (np.arange(ntimes)/float(ntimes)) * height
azimut_cam = (np.arange(ntimes)/float(ntimes)) * azimut

# Moving slice
if moving_h_slice:
	levels.insert(1, mhs_origin)
	mhs_delta = mhs_finish - mhs_origin
if moving_v_slice:
	mvs_delta_lon = mvs_finish_lon - mvs_origin_lon
	mvs_delta_lat = mvs_finish_lat - mvs_origin_lat

# Zooming parameters
dx_zoom = delta_lon * 0.1
dy_zoom = delta_lat * 0.1
dz_zoom = (z_lims[0] - z_lims[-1]) * 0.1

# Loop
t_roms = -1
frame = 0

for time_p in np.arange(t0, ntimes, tsep):

	
	t_roms_prev = t_roms
	t_roms = time_p // nsteps # Current step for the ROMS fields.

	t4int = time_p # Current step for the interpolated ROMS fields.
	
	## Read Fields
	read_now = t_roms != t_roms_prev
	if fields_calc:
		if read_now:
			with Dataset(roms_file) as nc:
				zeta[0] = nc.variables['zeta'][t_roms][region_eddy]
				zeta[1] = nc.variables['zeta'][t_roms + 1][region_eddy]
				if want['pt'] or want['Density']:
					pt[0] = nc.variables['temp'][t_roms][region3d_eddy]
					pt[1] = nc.variables['temp'][t_roms + 1][region3d_eddy]
				if want['salt'] or want['Density']:
					salt[0] = nc.variables['salt'][t_roms][region3d_eddy]
					salt[1] = nc.variables['salt'][t_roms + 1][region3d_eddy]
				if want['w'] or want['uv'] or want['h_vorticity']:
					#u[0] = nc.variables['u'][t_roms][rg3d_u_eddy]
					#v[0] = nc.variables['v'][t_roms][rg3d_v_eddy]
					u[0, :, :, :-1] = nc.variables['u'][t_roms][rg3d_u_eddy]
					v[0, :, :-1, :] = nc.variables['v'][t_roms][rg3d_v_eddy]	
					#u[1] = nc.variables['u'][t_roms + 1][rg3d_u_eddy]
					#v[1] = nc.variables['v'][t_roms + 1][rg3d_v_eddy]
					u[1, :, :, :-1] = nc.variables['u'][t_roms + 1][rg3d_u_eddy]
					v[1, :, :-1, :] = nc.variables['v'][t_roms + 1][rg3d_v_eddy]
				if want['w']: 	# This is not compatible with z_fields, by the moment.
					omega[0] = nc.variables['omega'][t_roms][region3d_eddy]
					omega[1] = nc.variables['omega'][t_roms + 1][region3d_eddy]


	#################################################################
	### Conversion from omega to w ###
	#################################################################

	# This may go after the conversion to z levels... No, I don't think so.

			if want['w'] and not z_fields: # omega_to_w is still not prepared for z_fields.
				print '  - Conversion from omega to w. '
				print
				w[0] = omega_to_w(omega[0], u[0, :, :, :-1], v[0, :, :-1, :], zeta[0], h_neg[region_eddy], bound=15.)[:]
				w[1] = omega_to_w(omega[1], u[1, :, :, :-1], v[1, :, :-1, :], zeta[1], h_neg[region_eddy], bound=15.)[:]

	#################################################################
	### Pressure and density ###
	#################################################################
			'''		
			print '  - Denisty.'

			if fade_initial or other_fields:

				# Pressure
				print '  - Pressure.'
				press[:] = gsw.p_from_z(np.repeat(np.asarray(lev_fie_chosen), y_eddy.size), 
					np.concatenate(tuple([y_eddy.flatten() for i in lev_fie_chosen]))).reshape(press.shape)[:]
			'''
			if not zrt_indep_t:
				zrt = octant.depths.get_zrho(Vtransform, Vstretching, km, theta_s, theta_b, 
					                    -h_neg[region_eddy], hc, zeta=zeta[0], Hscale=3)
			'''
			for t in np.arange(fieldtimes):
				dens[t, :] = gsw.rho(salt[t, :].flatten(), 
							gsw.CT_from_pt(salt_r[t, :].flatten(), pt[t, :].flatten()), 
							press[:].flatten()).reshape(dens.shape[1:])
				
				dens[t] = rho_eos_duko(pt, salt, zrt, rho0 = 1027.40898956694)

			fields_s.append(dens)
			print
			'''	
	#################################################################
	### Convert fields from sigma to z levels ###
	#################################################################

			"""
			Interoplate the fields in the vertical to obtain z_levels fields from the sigma_levels ones.
			The wanted number of z_levels is set up in the Control Pannel. In principle, values under topography may be either zero or masked.
			The z levels range within the limits of the levels list.
			"""

			print '  - Interpolating from sigma to z levels. '

			for i, field in enumerate(fields_s):
				if zrt_indep_t:
					for kz, z in enumerate(levels_field):
						for ks in np.arange(roms_levels - 1):

							if zrtm1[kz].size > 0:
								for t in np.arange(2):
									fields_z[i][t, kz][in_dz[kz]] = field[t, ks][in_dz[kz]] + (field[t, ks + 1][in_dz[kz]] - field[t, ks][in_dz[kz]]) * weigh[kz]
				else:
					for kz, z in enumerate(levels_field):
						for ks in np.arange(km - 1):

							zrtm1 = zrt[ks]
							zrtp1 = zrt[ks + 1]

							in_dz = np.where((zrtm1 < z) & (zrtp1 >= z))

							zrtm1 = zrtm1[in_dz]
							zrtp1 = zrtp1[in_dz]

							if zrtm1.size > 0:
								weigh = abs((z - zrtm1) / (zrtp1 - zrtm1))
								for t in np.arange(2):
									fields_z[i + 1][t, kz][in_dz] = field[t, ks][in_dz] + (field[t, ks + 1][in_dz] - field[t, ks][in_dz]) * weigh
									# i + 1 because fields_z has zeta in the first position and fields_s does not.
									if np.where(np.isnan(fields_z[i + 1][t, kz][in_dz]))[0].size > 1:
										print 'Here we have nans.', i, t, kz, in_dz
									if np.where(fields_z[i + 1][t, kz][in_dz] < 0.)[0].size > 1:
										print 'Here we have bad values.', i, t, kz, in_dz

	#############################################################
	### Create new grid and interpolate ALL the fields onto it
	#############################################################

			fields_rect_eddy = new_grid(son_eddy, sox_eddy, san_eddy, sax_eddy, 
						x_eddy, y_eddy, fields_z, new_shape=new_shape, only_xy=False)[2]

			# Assign values to the fields_r
			zeta_rect[:] = fields_rect_eddy[0][:]
			i = 1
			if want['w']:
				w_r[:] = fields_rect_eddy[i][:]
				i += 1
			if want['w'] or want['uv'] or want['h_vorticity']:
				u_r[:] = fields_rect_eddy[i][:]
				i += 1
				v_r[:] = fields_rect_eddy[i][:]
				i += 1
			if want['pt'] or want['Density']:
				pt_r[:] = fields_rect_eddy[i][:]
				if thereis_land:
					combined_mask_6[:] = np.ma.mask_or(pt_r.mask, np.ma.masked_where(pt_r < pt_z[np.where(pt_z != pt_z.min())].min(), pt_r).mask)
					pt_r[:] = np.ma.array(pt_r.data, mask=combined_mask_6[:])
				i += 1
			if want['salt'] or want['Density']:
				salt_r[:] = fields_rect_eddy[i][:]
				if thereis_land:
					combined_mask_6[:] = np.ma.mask_or(salt_r.mask, np.ma.masked_where(salt_r < pt_z[np.where(salt_z != salt_z.min())].min(), salt_r).mask)
					salt_r[:] = np.ma.array(salt_r.data, mask=combined_mask_6[:])
				i += 1
		
			if bathymetry:
				field_list_basin = [h_neg_basin]

				x_basin_old = x_basin[:]
				y_basin_old = y_basin[:]
				x_basin, y_basin, fields_rect_basin = new_grid(son_basin, sox_basin, san_basin, sax_basin, 
															x_basin, y_basin, field_list_basin, new_shape=(100, 100))
				h_neg_basin = fields_rect_basin[0]
		
	###############################################
	###	Calculate other fields. ###
	###############################################
			
			print '*** Calculating other fields. ***'
			print

			if want['h_vorticity']:
				# Vorticity
				print '  - Vertical component of the Vorticity.'
				dx = 1000. 	# This is not really true. Even less with the new_shape_f, but see if it works fine.
				for t in np.arange(u_r.shape[0]):
					for k in np.arange(u_r.shape[1]):
						hvort[t, k, :] = np.gradient(v_r[t, k, :], dx)[1] - np.gradient(u_r[t, k, :], dx)[0]
				badinds = np.where(hvort > 1e3) # Check for bad values
				hvort[badinds] = 0.

			
			if want['Density']:
				# Pressure
				print '  - Pressure.'
				press_r[:] = gsw.p_from_z(np.repeat(np.asarray(lev_fie_chosen), y_eddy_r.size), 
					np.concatenate(tuple([y_eddy_r.flatten() for i in lev_fie_chosen]))).reshape(press_r.shape)[:]
				# Density
				print '  - Density.'
				for t in np.arange(fieldtimes):
					dens_r[t] = gsw.rho(salt_r[t].flatten(), 
								gsw.CT_from_pt(salt_r[t].flatten(), pt_r[t].flatten()), 
								press_r[:].flatten()).reshape(dens_r.shape[1:])
					#dens[t] = rho_eos_duko(pt_r, salt_r)

				print

	###################################################################
	### Start what was in the old loop for making and saving frames ###
	###################################################################

	firstbar_fl = True

	# Timing variables:
	day = time_p // framesperday
	hour = time_p - 24 * day
	
	# Figure and plots:
	#figsize=(18, 11)
	fig = plt.figure(figsize=figsize)
	ax = Axes3D(fig)
	ax.set_xlabel('Longitude (degrees E)')
	ax.set_ylabel('Latitude (degrees N)')
	ax.set_zlabel('Depth (m)')

	if zoom:
		x_zoom = dx_zoom * 100 * (1 - time_p / ntimes)
		y_zoom = dy_zoom * 100 * (1 - time_p / ntimes)
		z_zoom = dz_zoom * 100 * (1 - time_p / ntimes)
		ax.set_xlim(x_lims[0] - x_zoom, x_lims[1] + x_zoom)
		ax.set_ylim(y_lims[0] - y_zoom, y_lims[1] + y_zoom)
		ax.set_zlim(levels[0] - z_zoom, levels[-1] + z_zoom)
	else:
		ax.set_xlim(x_lims)
		ax.set_ylim(y_lims)
		ax.set_zlim(z_lims)
	if show_time:
		#ax.set_title('Day %i, hour %i.'%(day, hour))
		ax.set_title('Day %i'%day)
	if move_cam:
		ax.view_init(height_cam[time_p], azimut_cam[time_p])
	else:
		ax.view_init(height, azimut)
	#ax.view_init(40, 30)
	if bathymetry: include_bathymetry(_lon, _lat)

	# Choose only the particles that are unmasked

	lons_now = lons_domain[time_p].compressed()
	lats_now = lats_domain[time_p].compressed()
	depths_now = depths_domain[time_p].compressed()

	
	if fields_calc:

		#print 'Interpolation of the fields' + cont_field + 'and' + fill_field + 'in time.'

		t4f = t_roms # Current step for the ROMS fields.
		t4int = time_p # Current step for the interpolated ROMS fields.

		## Interpolate the fields in time

		if show_black_field:
			if cont_field == 'salt':
				f2plot1[:] = salt_r[:]
				mean1 = salt_mean
				max1 = salt_maxim
				min1 = salt_minim
			if cont_field == 'pt':
				f2plot1[:] = pt_r[:]
				mean1 = pt_mean
				max1 = pt_maxim
				min1 = pt_minim

		if show_colored_field:
			if fill_field == 'salt':
				f2plot2[:] = salt_r[:]
				mean2 = salt_mean
				max2 = salt_maxim
				min2 = salt_minim
			if fill_field == 'pt':
				f2plot2[:] = pt_r[:]
				mean2 = pt_mean
				max2 = pt_maxim
				min2 = pt_minim

		f2plot1_int[:] = interp_field_1D(f2plot1[:], nnew, which=time_p-nisteps*t4f)[:]
		f2plot2_int[:] = interp_field_1D(f2plot2[:], nnew, which=time_p-nisteps*t4f)[:]

		if want['w']:
			w_z_int = interp_field_1D(w_z[t4f:t4f+2,:], nnew, which=time_p - nisteps*t4f)
			
		#if time_p == 1: sys.exit()
		# Moving slice

		#print 'Moving slices. Horizontal | Vertical:', moving_h_slice, moving_v_slice

		if moving_h_slice:
			moving_level = mhs_origin + abs(np.sin((time_p/ntimes) * 1. * np.pi)) * mhs_delta
			moving_level = mhs_origin + slice_movement(time_p, mov_type) * mhs_delta
			levels[1] = moving_level
		if moving_v_slice:
			if which_slices[0] and not which_slices[1]:
				wall_offset_tmp = mvs_origin_lon + slice_movement(time_p, mov_type) * mvs_delta_lon
				#jj = int(x_eddy_r.shape[1] * (1 - ((wall_offset_tmp - son_eddy) / delta_lon)))
				#jj = int(x_eddy_r.shape[1] * ((wall_offset_tmp - son_eddy) / delta_lon))
				jj = int(x_eddy_r.shape[1] * ((wall_offset_tmp - son_eddy) / delta_lon))
				wo_j = 1 if jj == 0 else jj
				wall_offset = x_eddy_r[0, wo_j - 1]
				#print time_p, (wall_offset_tmp - son_eddy) / delta_lon, jj, wo_j
			elif which_slices[1] and not which_slices[0]:
				wall_offset_tmp = mvs_origin_lat + slice_movement(time_p, mov_type) * mvs_delta_lat
				#ii = int(y_eddy_r.shape[0] * (1 - ((wall_offset_tmp - san_eddy) / delta_lat)))
				#ii = int(y_eddy_r.shape[0] * ((wall_offset_tmp - san_eddy) / delta_lat))
				ii = int(y_eddy_r.shape[0] * ((wall_offset_tmp - san_eddy) / delta_lat))
				wo_i = 1 if ii == 0 else ii
				wall_offset = y_eddy_r[wo_i - 1, 0]
			elif which_slices[0] and which_slices[1]:
				t4lon_slice = ntimes / 2
				if time_p < t4lon_slice:
					wall_offset_tmp_lon = mvs_origin_lon + slice_movement(time_p, mov_type) * mvs_delta_lon
					#jj = int(x_eddy_r.shape[1] * (1 - ((wall_offset_tmp_lon - sox_eddy) / delta_lon)))
					#jj = int(x_eddy_r.shape[1] * ((wall_offset_tmp_lon - sox_eddy) / delta_lon))
					jj = int(x_eddy_r.shape[1] * ((wall_offset_tmp_lon - sox_eddy) / delta_lon))
					wo_j = 1 if jj == 0 else jj
					wall_offset = x_eddy_r[0, wo_j - 1]
				elif time_p >= t4lon_slice:
					wall_offset_tmp_lat = mvs_origin_lat + slice_movement(time_p, mov_type) * mvs_delta_lat
					#ii = int(y_eddy_r.shape[0] * (1 - ((wall_offset_tmp_lat - sax_eddy) / delta_lat)))
					#ii = int(y_eddy_r.shape[0] * ((wall_offset_tmp_lat - sax_eddy) / delta_lat))
					ii = int(y_eddy_r.shape[0] * ((wall_offset_tmp_lat - sax_eddy) / delta_lat))
					wo_i = 1 if ii == 0 else ii
					wall_offset = y_eddy_r[wo_i - 1, 0]


		if fields_plot:

			#print 'Plotting the fields'

			firstbar_fie = True

			# Separate the particles with the slice as a separatrix
			
			if cube:
				sep_lev = levels[0]
			else:
				sep_lev = levels[1]

			if horiz_slices or cube:
				mask_above = np.ma.masked_less(depths_now, sep_lev).mask
				mask_below = np.ma.masked_greater(depths_now, sep_lev).mask
			elif vertic_slices:
				if which_slices[0]:
					mask_above = np.ma.masked_less(lons_now, wall_offset).mask
					mask_below = np.ma.masked_greater(lons_now, wall_offset).mask
				elif which_slices[1]:
					mask_above = np.ma.masked_less(lats_now, wall_offset).mask
					mask_below = np.ma.masked_greater(lats_now, wall_offset).mask
			else:
				depths_above = depths_now
				depths_below = depths_now
				lons_above = lons_now
				lons_below = lons_now
				lats_above = lats_now
				lats_below = lats_now
			# Check this
			depths_above = np.ma.array(depths_now, mask=mask_above).compressed()
			depths_below = np.ma.array(depths_now, mask=mask_below).compressed()
			lons_above = np.ma.array(lons_now, mask=mask_above).compressed()
			lons_below = np.ma.array(lons_now, mask=mask_below).compressed()
			lats_above = np.ma.array(lats_now, mask=mask_above).compressed()
			lats_below = np.ma.array(lats_now, mask=mask_below).compressed()
			
			# Indices for the levels that are already in levels_field
			old_levels_inds = []
			new_levels_inds = []
			for i, level in enumerate(levels):
				if level in levels_field:
					old_levels_inds.append(i)
				else:
					new_levels_inds.append(i)

			if fade_slices:
				alpha_slices = 1 - time_p / ntimes

			if horiz_slices:
				create_slices(f2plot1_int, f2plot2_int)

			for i, level in enumerate(levels):	
				# Scatter the particles below (behind) the slice first
				if show_particles and time_p != 0:
					if fade_method is not None:
						trackscatter_all(lons_below, lats_below, depths_below, x_eddy, y_eddy, firstbar_fl, colorby=colorby,
										 fade_field=fade_field, fade_method=fade_method, what2do='plot', zorder=1)
					else:
						trackscatter_all(lons_below, lats_below, depths_below, x_eddy, y_eddy, firstbar_fl, colorby=colorby,
										 color='Green', fade_method=fade_method, what2do='plot', zorder=1)
				
				# Plot the slices. Avoid limits of the level list.
				if i != 0 and level != levels[-1]:

					if horiz_slices:
						draw_fields(field_slice[i, :], field_slice2[i, :], v_raw=[f2plot1_int, f2plot2_int], what2draw='Slices')


					if cube:
						if moving_v_slice:
							if which_slices[0] and not which_slices[1]:
								draw_fields(f2plot1_int[:, :, :wo_j], f2plot2_int[:, :, :wo_j], x_eddy_r[0, :wo_j], y_eddy_r[:,0], what2draw='Cube',
									offsets=(lfc[0], y_eddy_r.max(), wall_offset))

							elif which_slices[1] and not which_slices[0]:
								draw_fields(f2plot1_int[:, :wo_i], f2plot2_int[:, :wo_i], x_eddy_r[0,:], y_eddy_r[:wo_i, 0], what2draw='Cube',
									offsets=(lfc[0], wall_offset, x_eddy_r.max()))

							elif which_slices[0] and which_slices[1]:
								if time_p < t4lon_slice:
									draw_fields(f2plot1_int[:, :, :wo_j], f2plot2_int[:, :, :wo_j], x_eddy_r[0, :wo_j], y_eddy_r[:,0], what2draw='Cube',
										offsets=(lfc[0], y_eddy_r.max(), wall_offset))

								elif time_p >= t4lon_slice:
									draw_fields(f2plot1_int[:, :wo_i], f2plot2_int[:, :wo_i], x_eddy_r[0,:], y_eddy_r[:wo_i, 0], what2draw='Cube',
										offsets=(lfc[0], wall_offset, x_eddy_r.max()))

						else:
							draw_fields(f2plot1_int, f2plot2_int, x_eddy_r[0,:], y_eddy_r[:,0], what2draw='Cube',
								offsets=(lfc[0], y_eddy_r.max(), x_eddy_r.max()))


					if vertic_slices:
						wall_slices(which=which_slices)

					firstbar_fie = False
				
				# Scatter the particles above (situated before) the slice
				if show_particles and time_p != 0:
					firstbar_fl = False
					if fade_method is not None:
						trackscatter_all(lons_above, lats_above, depths_above, x_eddy, y_eddy, firstbar_fl, colorby=colorby,
										 fade_field=fade_field, fade_method=fade_method, what2do='scatter', zorder=3)
					else:
						trackscatter_all(lons_above, lats_above, depths_above, x_eddy, y_eddy, firstbar_fl, colorby=colorby,
										 color='Red', fade_method=fade_method, what2do='scatter', zorder=3)

		#t4f_old = t4f

	nlats_now = lats_domain[time_p].count()
	nlons_now = lons_domain[time_p].count()

	# Scatter the particles alone
	if show_particles and not fields_plot:
		if fade_method is not None:
			trackscatter_all(lons_now, lats_now, depths_now, x_eddy, y_eddy, firstbar_fl, colorby=colorby,
							 fade_field=fade_field, fade_method=fade_method)
		else:
			trackscatter_all(lons_now, lats_now, depths_now, x_eddy, y_eddy, firstbar_fl, colorby=colorby,
							 fade_method=fade_method)

	#plt.tight_layout()
	#plt.gca().invert_zaxis()
	#plt.invert_zaxis()
	ax.invert_zaxis()	
	if not show_axes: plt.axis('off')
	image_name = folder_path + anim_name + '_Frame_%04i'%frame + '.png'
	fig.savefig(image_name, dpi=dpi)
	#plt.show()
	plt.close(fig)

	## Resizing a large image to the wanted resolution
	size = 1920, 1080
	image = Image.open(image_name)
	im_resized = image.resize(size, Image.ANTIALIAS)
	im_resized.save(image_name, "PNG")
	print 'Frame created.', anim_nu, frame, lats_now.size #, nlons_now
	frame += 1

	#print '-------------------------------------------------------------'
