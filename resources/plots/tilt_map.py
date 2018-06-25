#!/usr/bin/env python

from __future__ import print_function

import math
import numpy

from os.path import expandvars

import matplotlib
matplotlib.use("PDF")

fig_size = [11.7,8.3] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 12,
        'text.fontsize': 12,
        'legend.fontsize': 8,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': True,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

import pylab
import scipy
import scipy.interpolate


from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units


def addAnnotationToPlot(plot, text, loc=1, size=6.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size), frameon=True,
                      loc=loc,
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)

wlens=numpy.linspace(300.,600.,num=100)



detectorCenterDepth = 1948.07*I3Units.m


tilt_distance_from_origin_in_tilt_dir = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/ice/TILT_data/tilt.par"), unpack=True)[1]*I3Units.m

tilt_dat = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/ice/TILT_data/tilt.dat"), unpack=True)
tilt_zcoords = (detectorCenterDepth-tilt_dat[0])[::-1]

tilt_zcoord_diffs = tilt_zcoords[1:]-tilt_zcoords[:-1]
z_offset = numpy.mean(tilt_zcoord_diffs)

tilt_shift = []
for i in range(len(tilt_distance_from_origin_in_tilt_dir)):
    tilt_shift.append(tilt_dat[i+1][::-1])
tilt_shift = numpy.array(tilt_shift)

# print "dist:", tilt_distance_from_origin_in_tilt_dir
# print "zcoords:", tilt_zcoords
# #print "zcoord_diffs:", tilt_zcoord_diffs
# print "shift:", tilt_shift

# direction of ice tilt
tiltAngle = 225.*I3Units.deg

# (x,y) position of string 50
lnx = numpy.cos(tiltAngle)
lny = numpy.sin(tiltAngle)


def zshift(rx, ry, rz):
    #print "eval at", rx, ry, rz

    z = (rz-tilt_zcoords[0])/z_offset
    k = min(max(numpy.floor(z), 0), len(tilt_zcoords)-2)
    l=k+1

    nr = lnx*rx + lny*ry

    for j in range(1,len(tilt_distance_from_origin_in_tilt_dir)):
        if (nr < tilt_distance_from_origin_in_tilt_dir[j]) or (j==len(tilt_distance_from_origin_in_tilt_dir)-1):
            frac_at_lower = (tilt_distance_from_origin_in_tilt_dir[j] - nr  )/(tilt_distance_from_origin_in_tilt_dir[j] - tilt_distance_from_origin_in_tilt_dir[j-1])
            frac_at_upper = (nr - tilt_distance_from_origin_in_tilt_dir[j-1])/(tilt_distance_from_origin_in_tilt_dir[j] - tilt_distance_from_origin_in_tilt_dir[j-1])

            val_at_lower = (tilt_shift[j-1][l]*(z-k) + tilt_shift[j-1][k]*(l-z))
            val_at_upper = (tilt_shift[j]  [l]*(z-k) + tilt_shift[j]  [k]*(l-z));

            return ( val_at_upper * frac_at_upper +
                     val_at_lower * frac_at_lower )

zshift_vectorized = numpy.vectorize(zshift)


iceTiltCLSim = clsim.util.GetIceTiltZShift()
zshiftCLSim_vectorized = numpy.vectorize(lambda x,y,z: iceTiltCLSim.GetValue(x,y,z))


def plotAtDepth(plot, scanDepth):
    scanZ = detectorCenterDepth-scanDepth

    delta = 10.
    x = numpy.arange(-600.0, 600.0+delta, delta)
    y = numpy.arange(-600.0, 600.0+delta, delta)
    X, Y = numpy.meshgrid(x, y)

    Z_reference = zshift_vectorized(X, Y, scanZ)
    Z = zshiftCLSim_vectorized(X, Y, scanZ)

    comp = numpy.array([numpy.ravel(Z), numpy.ravel(Z_reference)]).T
    for a, b in comp:
        diff = numpy.abs(a-b)
        if diff > 1e-10:
            raise RuntimeError("reference and c++ results differ!")
    del comp

    level_delta = 2.
    level_from = -30.
    level_to = 50.

    levels = numpy.arange(level_from,level_to+level_delta,level_delta)

    extent = [x[0], x[-1], y[0], y[-1]]
    im = plot.imshow(
        Z.T,
        extent=extent,
        interpolation='nearest',
        origin='lower',
        cmap=matplotlib.pyplot.cm.hot,
        norm=matplotlib.colors.Normalize(vmin=level_from, vmax=level_to))
    # cbar = fig.colorbar(im, ax=plot)

    CS = plot.contour(X, Y, Z, levels=levels, colors='k')
    plot.clabel(CS, inline=1, fontsize=10, fmt=r'$%+1.1f\,\mathrm{m}$')


    plot.grid(True)

    plot.set_xlabel("x $[\\mathrm{m}]$")
    plot.set_ylabel("y $[\\mathrm{m}]$")
    addAnnotationToPlot(plot, r"depth=$%1.0f\,\mathrm{m}$ / z=$%1.0f\,\mathrm{m}$" % (scanDepth, scanZ), loc=2)

    plot.set_xlim(-600.,600.)
    plot.set_ylim(-600.,600.)


fig = pylab.figure()
fig.subplots_adjust(left=0.06, bottom=0.055, top=0.96, right=0.98, hspace=0.15, wspace=0.25)

ax = fig.add_subplot(2, 3, 1)
bx = fig.add_subplot(2, 3, 3)
cx = fig.add_subplot(2, 3, 2)
dx = fig.add_subplot(2, 3, 4)
ex = fig.add_subplot(2, 3, 5)
fx = fig.add_subplot(2, 3, 6)

plotAtDepth(bx, 1250.*I3Units.m)
plotAtDepth(cx, 1550.*I3Units.m)
plotAtDepth(dx, 1850.*I3Units.m)
plotAtDepth(ex, 2150.*I3Units.m)
plotAtDepth(fx, 2400.*I3Units.m)

outfileName = "tilt_map.pdf"
pylab.savefig(outfileName, transparent=False)
print("wrote", outfileName)

