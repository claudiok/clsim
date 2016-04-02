#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use("PDF")


from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)

import os

options,args = parser.parse_args()

if len(args) != 0:
    parser.error("wrong number of options")


outfile = "timing_distributions_SpiceLea_tiltOnOff_anisotropyOnOff.pdf"

fig_size = [11.7,8.3] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 6,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

def addAnnotationToPlot(plot, text, loc=1, size=8., rotation=0.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size, rotation=rotation),
                      frameon=True,
                      loc=loc, # 1=='upper right'
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)

def addAnnotationToPlot(plot, text, loc=1, size=6.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size), frameon=True,
                      loc=loc,
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)

import math
import numpy
import pylab
import scipy
import scipy.interpolate
import scipy.integrate

import tables
import numpy

from icecube import icetray, dataclasses, clsim
from I3Tray import I3Units

def getTimesAndPositions(filename, DOMs):
    xPosForDOM = []
    yPosForDOM = []
    zPosForDOM = []

    timesForDOM = []

    h5file = tables.openFile(filename=filename, mode='r')
    allTimes = h5file.root.MCPESeriesMap.cols.time[:]
    allOMs = h5file.root.MCPESeriesMap.cols.om[:]
    allStrings = h5file.root.MCPESeriesMap.cols.string[:]
    allX = h5file.root.MCPESeriesMap.cols.x[:]
    allY = h5file.root.MCPESeriesMap.cols.y[:]
    allZ = h5file.root.MCPESeriesMap.cols.z[:]

    for string, dom in DOMs:
        timesForDOM.append(allTimes[(allOMs==dom) & (allStrings==string)])

        xPosForThisDOM = numpy.unique(allX[(allOMs==dom) & (allStrings==string)])
        yPosForThisDOM = numpy.unique(allY[(allOMs==dom) & (allStrings==string)])
        zPosForThisDOM = numpy.unique(allZ[(allOMs==dom) & (allStrings==string)])

        if (len(xPosForThisDOM) > 1) or (len(yPosForThisDOM) > 1) or (len(zPosForThisDOM) > 1):
            raise RuntimeError("DOM positions are not unique!")

        xPosForDOM.append(xPosForThisDOM[0])
        yPosForDOM.append(yPosForThisDOM[0])
        zPosForDOM.append(zPosForThisDOM[0])

    numEvents = len(h5file.root.__I3Index__.MCPESeriesMap.cols.exists[:])

    emitterPosX = numpy.unique(h5file.root.MCMostEnergeticInIce.cols.x[:])
    emitterPosY = numpy.unique(h5file.root.MCMostEnergeticInIce.cols.y[:])
    emitterPosZ = numpy.unique(h5file.root.MCMostEnergeticInIce.cols.z[:])

    if len(emitterPosX) != 1: raise RuntimeError("all emitters/particles need to be identical!")
    if len(emitterPosY) != 1: raise RuntimeError("all emitters/particles need to be identical!")
    if len(emitterPosZ) != 1: raise RuntimeError("all emitters/particles need to be identical!")

    emitterPos = numpy.array([emitterPosX[0], emitterPosY[0], emitterPosZ[0]])

    h5file.close()

    return (timesForDOM, xPosForDOM, yPosForDOM, zPosForDOM, numEvents, emitterPos)


OMKeys = [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (8, 2)]

filenames = [
    "test_events_clsim_lea.i3.hdf5",
    "test_events_clsim_lea_notilt.i3.hdf5",
    "test_events_clsim_lea_noanisotropy.i3.hdf5",
    "test_events_clsim_lea_notilt_noanisotropy.i3.hdf5",

    "test_events_ppc_lea.i3.hdf5",
    "test_events_ppc_lea_notilt.i3.hdf5",
    "test_events_ppc_lea_noanisotropy.i3.hdf5",
    "test_events_ppc_lea_notilt_noanisotropy.i3.hdf5",
]
colors = [
    'k',
    'r',
    'g',
    'b',
    'k',
    'r',
    'g',
    'b',
    ]
linestyles = [
    '-',
    '-',
    '-',
    '-',
    '--',
    '--',
    '--',
    '--',
]
labels = [
    'clsim std',
    'clsim no tilt',
    'clsim no aniso.',
    'clsim no tilt/no aniso.',
    'PPC std',
    'PPC no tilt',
    'PPC no aniso.',
    'PPC no tilt/no aniso.',
]
show = [
    True,     # clsim std
    True,     # clsim no tilt
    True,     # clsim no aniso.
    True,     # clsim no tilt/no aniso.

    True,     # PPC std
    True,     # PPC no tilt
    True,     # PPC no aniso.
    True,     # PPC no tilt/no aniso.
]


print("loading data..")

DOMpositionsX = numpy.ones(len(OMKeys)) * float('NaN')
DOMpositionsY = numpy.ones(len(OMKeys)) * float('NaN')
DOMpositionsZ = numpy.ones(len(OMKeys)) * float('NaN')

timesForFilename = []
numEventsForFilename = []

emitterPos = None

for filename in filenames:
    print("reading", filename)

    times, xPos, yPos, zPos, numEvents, thisEmitterPos = getTimesAndPositions(filename, OMKeys)

    if emitterPos is None:
        emitterPos = thisEmitterPos
    else:
        if thisEmitterPos[0] != emitterPos[0] or thisEmitterPos[1] != emitterPos[1] or thisEmitterPos[2] != emitterPos[2]:
            raise RuntimeError("input files cannot have emitting particles in different positions!")

    timesForFilename.append(times)
    numEventsForFilename.append(numEvents)

    for i in range(len(OMKeys)):
        key = OMKeys[i]
        if (numpy.isnan(DOMpositionsX[i])):
            DOMpositionsX[i] = xPos[i]
        else:
            if DOMpositionsX[i] != xPos[i]:
                print("got:", xPos)
                print("expected:", DOMpositionsX)

                raise RuntimeError("files have inconsistent DOM positions (x)")
        if (numpy.isnan(DOMpositionsY[i])):
            DOMpositionsY[i] = yPos[i]
        else:
            if DOMpositionsY[i] != yPos[i]:
                print("got:", xPos)
                print("expected:", DOMpositionsX)

                raise RuntimeError("files have inconsistent DOM positions (y)")
        if (numpy.isnan(DOMpositionsZ[i])):
            DOMpositionsZ[i] = zPos[i]
        else:
            if DOMpositionsZ[i] != zPos[i]:
                print("got:", xPos)
                print("expected:", DOMpositionsX)

                raise RuntimeError("files have inconsistent DOM positions (z)")



print("done.")



####
print("plotting..")


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

subplots = [
    fig.add_subplot(3, 3, 2),
    fig.add_subplot(3, 3, 3),
    fig.add_subplot(3, 3, 6),
    fig.add_subplot(3, 3, 9),
    fig.add_subplot(3, 3, 8),
    fig.add_subplot(3, 3, 7),
    fig.add_subplot(3, 3, 4),
    fig.add_subplot(3, 3, 1),
]


def plotHistogram(plot, times, weights=None, color='k', linestyle='-', label=None):
    the_range=(500.,2500.)
    num_bins=200

    hist, bin_edges = numpy.histogram(times, weights=weights, bins=num_bins, range=the_range)
    plot.semilogy( (bin_edges[1:]+bin_edges[:-1])/2., hist, color=color, linestyle=linestyle, label=label)

for i in range(len(timesForFilename)):
    if not show[i]: continue

    timesForDOM = timesForFilename[i]
    filename = filenames[i]
    label = labels[i]
    linestyle = linestyles[i]
    color = colors[i]
    numEventsInFile = numEventsForFilename[i]

    for j, times in enumerate(timesForDOM):
        subplot = subplots[j]
        weights = numpy.ones(len(times)) / float(numEventsInFile)

        plotHistogram(subplot, times, weights=weights, color=color, linestyle=linestyle, label=label)


for subplot in subplots:
    subplot.grid(True)
    subplot.set_xlim(500., 2500.)
    subplot.set_ylim(3e-1, 3e1)
    subplot.legend(loc='upper right')

    subplot.set_xlabel(r"$t_\mathrm{hit;MC}$ [$\mathrm{ns}$]")
    subplot.set_ylabel(r"$N_\mathrm{hit;MC}$")


centerPlotRangeX = [-160.+emitterPos[0],160.+emitterPos[0]]
centerPlotRangeY = [-160.+emitterPos[1],160.+emitterPos[1]]



ax = fig.add_subplot(3, 3, 5)
ax.set_aspect('equal')



# plot a tilt map for reference
detectorCenterDepth = 1948.07*I3Units.m

iceTiltCLSim = clsim.util.GetIceTiltZShift()
zshiftCLSim_vectorized = numpy.vectorize(lambda x,y,z: iceTiltCLSim.GetValue(x,y,z))

scanZ = emitterPos[2]
scanDepth = detectorCenterDepth-scanZ

delta = 1.
x = numpy.arange(centerPlotRangeX[0], centerPlotRangeX[1]+delta, delta)
y = numpy.arange(centerPlotRangeY[0], centerPlotRangeY[1]+delta, delta)
X, Y = numpy.meshgrid(x, y)
Z = zshiftCLSim_vectorized(X, Y, scanZ)

level_delta = 4.
level_from = -30.
level_to = 50.

levels = numpy.arange(level_from,level_to+level_delta,level_delta)

extent = [x[0], x[-1], y[0], y[-1]]
im = ax.imshow(
    Z.T,
    extent=extent,
    interpolation='nearest',
    origin='lower',
    cmap=matplotlib.pyplot.cm.gray,
    norm=matplotlib.colors.Normalize(vmin=level_from, vmax=level_to))
# cbar = fig.colorbar(im, ax=plot)

CS = ax.contour(X, Y, Z, levels=levels, colors='k')
ax.clabel(CS, inline=1, fontsize=8, fmt=r'$%+1.1f\,\mathrm{m}$')

addAnnotationToPlot(ax, r"z=$%1.0f\,\mathrm{m}$ (depth=$%1.0f\,\mathrm{m}$)" % (scanZ, scanDepth), loc=2)






# add the middle plot
tilt_azimuth = 225. *numpy.pi/180.
tilt_dir_x = numpy.cos(tilt_azimuth)
tilt_dir_y = numpy.sin(tilt_azimuth)

# 216deg is the direction of tilt. we want the direction of flow here
anis_azimuth = (216. -90.) *numpy.pi/180.
anis_dir_x = numpy.cos(anis_azimuth)
anis_dir_y = numpy.sin(anis_azimuth)

arrow_magnitude=100.

ax.scatter(DOMpositionsX, DOMpositionsY, s=20, c='k', marker='o',zorder=2)
ax.scatter([emitterPos[0]], [emitterPos[1]], s=120, c='r', marker='o',zorder=2)
ax.scatter([emitterPos[0]], [emitterPos[1]], s=120/numpy.sqrt(2.)*0.85, c='k', marker='x', zorder=2)

c = matplotlib.patches.FancyArrowPatch(
    posA = (-tilt_dir_x*arrow_magnitude+emitterPos[0], -tilt_dir_y*arrow_magnitude+emitterPos[1]),
    posB = ( tilt_dir_x*arrow_magnitude+emitterPos[0],  tilt_dir_y*arrow_magnitude+emitterPos[1]),
    arrowstyle="simple",
    mutation_scale=100.,
    edgecolor='0.5',
    facecolor='0.5',
    alpha=0.5,
    zorder=1,
    )
ax.add_patch(c)


d = matplotlib.patches.FancyArrowPatch(
    posA = (-anis_dir_x*arrow_magnitude+emitterPos[0], -anis_dir_y*arrow_magnitude+emitterPos[1]),
    posB = ( anis_dir_x*arrow_magnitude+emitterPos[0],  anis_dir_y*arrow_magnitude+emitterPos[1]),
    arrowstyle="simple",
    mutation_scale=100.,
    edgecolor='b',
    facecolor='b',
    alpha=0.5,
    zorder=1,
    )
ax.add_patch(d)


ax.grid(True)
ax.set_xlim(centerPlotRangeX[0], centerPlotRangeX[1])
ax.set_ylim(centerPlotRangeY[0], centerPlotRangeY[1])
ax.set_xlabel(r"x [$\mathrm{m}$]")
ax.set_ylabel(r"y [$\mathrm{m}$]")

print("saving as {0}".format(outfile))
pylab.savefig(outfile, transparent=False)



