#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use("PDF")


from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
#parser.add_option("-f","--format",dest="format",help="format to output [hdf5, root, or csv]",default='hdf5')
#parser.add_option("-z","--compress",dest="compression",help="compression level",default=1,type=int)
#parser.add_option("-n","--frames",dest="nframes",help="number of frames to process",default=None,type=int)

import os

options,args = parser.parse_args()

if len(args) != 1:
    parser.error("You must supply an input filename")

infile = args[0]
outfile = os.path.basename(args[0]) + '.pdf'





fig_size = [11.7,8.3] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 6,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': False,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
# matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

def addAnnotationToPlot(plot, text, loc=1, size=8., rotation=0.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size, rotation=rotation),
                      frameon=True,
                      loc=loc, # 1=='upper right'
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

h5file = tables.openFile(filename=infile, mode='r')



def calcTimeResiduals(hitSeriesMap,
                      hitSeriesIndex,
                      muons,
                      muonIndex):

    if len(muonIndex) != len(hitSeriesIndex):
        raise RuntimeError("muon and hit indices have different lengths")

    muon_x = muons.x[:]
    muon_y = muons.y[:]
    muon_z = muons.z[:]
    muon_t = muons.time[:]
    muon_zenith = muons.zenith[:]
    muon_azimuth = muons.azimuth[:]

    muon_theta = math.pi - muon_zenith
    muon_phi = muon_azimuth - math.pi
    muon_rho = numpy.sin(muon_theta)
    muon_dx = muon_rho*numpy.cos(muon_phi)
    muon_dy = muon_rho*numpy.sin(muon_phi)
    muon_dz = numpy.cos(muon_theta)
    del muon_theta, muon_phi, muon_rho
    

    n_group = 1.35634
    n_phase = 1.3195
    theta_c = numpy.arccos(1./n_phase)
    tan_theta_c = numpy.tan(theta_c)
    sin_theta_c = numpy.sin(theta_c)
    c_light = 0.299792458 # m/ns
    c_photon = c_light/n_group

    hits_x = hitSeriesMap.x[:]
    hits_y = hitSeriesMap.y[:]
    hits_z = hitSeriesMap.z[:]
    hits_t = hitSeriesMap.time[:]
    hitsIndex_start = hitSeriesIndex.start[:]
    hitsIndex_stop = hitSeriesIndex.stop[:]

    muonIndex_start = muonIndex.start[:]
    muonIndex_stop = muonIndex.stop[:]

    time_residuals = numpy.zeros(len(hits_x))
    track_dca = numpy.zeros(len(hits_x))

    print("calculcating time residuals..")
    for i in range(len(muonIndex_start)):
        # print "i=", i
        # print "hitsIndex_start[i]=", hitsIndex_start[i]
        # print "hitsIndex_stop[i]=", hitsIndex_stop[i]

        x = hits_x[hitsIndex_start[i]:hitsIndex_stop[i]]
        y = hits_y[hitsIndex_start[i]:hitsIndex_stop[i]]
        z = hits_z[hitsIndex_start[i]:hitsIndex_stop[i]]
        t = hits_t[hitsIndex_start[i]:hitsIndex_stop[i]]
    
        m_x = muon_x[muonIndex_start[i]]
        m_y = muon_y[muonIndex_start[i]]
        m_z = muon_z[muonIndex_start[i]]
        m_dx = muon_dx[muonIndex_start[i]]
        m_dy = muon_dy[muonIndex_start[i]]
        m_dz = muon_dz[muonIndex_start[i]]
        m_t = muon_t[muonIndex_start[i]]
    
        v_x = x-m_x
        v_y = y-m_y
        v_z = z-m_z
    
        v_length = numpy.sqrt(v_x**2 + v_y**2 + v_z**2)
        dist_along = v_x*m_dx + v_y*m_dy + v_z*m_dz
        dist_to_track = numpy.sqrt(v_length**2 - dist_along**2)

        dist_to_emission = dist_along - dist_to_track/tan_theta_c
        cherenkov_dist = dist_to_track/sin_theta_c
    
        expected_time = m_t + dist_to_emission/c_light + cherenkov_dist/c_photon
    
        time_delay = t-expected_time

        time_residuals[hitsIndex_start[i]:hitsIndex_stop[i]] = time_delay
        track_dca[hitsIndex_start[i]:hitsIndex_stop[i]] = dist_to_track

    return (time_residuals, track_dca)

print("ppc...")
time_residuals_ppc, dca_ppc = calcTimeResiduals(h5file.root.MCHitSeriesMap.cols,
                                       h5file.root.__I3Index__.MCHitSeriesMap.cols,
                                       h5file.root.MCMostEnergeticMuon.cols,
                                       h5file.root.__I3Index__.MCMostEnergeticMuon.cols)

print("clsim..")
time_residuals_clsim, dca_clsim = calcTimeResiduals(h5file.root.MCPESeriesMap_clsim.cols,
                                         h5file.root.__I3Index__.MCPESeriesMap_clsim.cols,
                                         h5file.root.MCMostEnergeticMuon.cols,
                                         h5file.root.__I3Index__.MCMostEnergeticMuon.cols)


print("some more work..")
#
hits_string_ppc = h5file.root.MCHitSeriesMap.cols.string[:]
hits_string_clsim = h5file.root.MCPESeriesMap_clsim.cols.string[:]

hits_om_ppc = h5file.root.MCHitSeriesMap.cols.om[:]
hits_om_clsim = h5file.root.MCPESeriesMap_clsim.cols.om[:]


bincounts_ppc = numpy.bincount(hits_om_ppc[:][(dca_ppc>20.) & (hits_string_ppc==21)])
bincounts_clsim = numpy.bincount(hits_om_clsim[(dca_clsim>20.) & (hits_string_clsim==21)])
dom_numbers = list(range(0,len(bincounts_ppc)))


bincounts_ppc_string = numpy.bincount(hits_string_ppc[:][(dca_ppc>20.)])
bincounts_clsim_string = numpy.bincount(hits_string_clsim[(dca_clsim>20.)])
string_numbers = list(range(0,len(bincounts_ppc_string)))


####
print("plotting..")


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(2, 2, 1)
bx = fig.add_subplot(2, 2, 2)
cx = fig.add_subplot(2, 2, 3)
dx = fig.add_subplot(4, 2, 6)
ex = fig.add_subplot(4, 2, 8)

ax.scatter(dom_numbers, bincounts_ppc, marker='x', color='r', label='ppc')
ax.errorbar(dom_numbers, bincounts_ppc, yerr=numpy.sqrt(bincounts_ppc), xerr=0.5, fmt=None, ecolor='r')

ax.scatter(dom_numbers, bincounts_clsim, marker='x', color='b', label='clsim')
ax.errorbar(dom_numbers, bincounts_clsim, yerr=numpy.sqrt(bincounts_clsim), xerr=0.5, fmt=None, ecolor='b')

cx.scatter(numpy.array(dom_numbers, float), numpy.array(bincounts_clsim, float)/numpy.array(bincounts_ppc, float))




ax.set_xlim(-0.5,60.5)
ax.legend(loc='upper right')
ax.grid(True)
ax.set_xlabel("DOM number")
ax.set_ylabel("number of hits")

cx.set_xlim(-0.5,60.5)
cx.grid(True)
cx.set_xlabel("DOM number")



dx.scatter(string_numbers, bincounts_ppc_string, marker='x', color='r', label='ppc')
dx.errorbar(string_numbers, bincounts_ppc_string, yerr=numpy.sqrt(bincounts_ppc_string), xerr=0.5, fmt=None, ecolor='r')

dx.scatter(string_numbers, bincounts_clsim_string, marker='x', color='b', label='clsim')
dx.errorbar(string_numbers, bincounts_clsim_string, yerr=numpy.sqrt(bincounts_clsim_string), xerr=0.5, fmt=None, ecolor='b')

dx.set_xlim(-0.5,86.5)
dx.grid(True)
dx.set_xlabel("string number")


ex.scatter(numpy.array(string_numbers, float), numpy.array(bincounts_clsim_string, float)/numpy.array(bincounts_ppc_string, float))
ex.set_xlim(-0.5,86.5)
ex.grid(True)
ex.set_xlabel("string number")



#the_range=(-20.,20.)
the_range=(-10.,50.)
num_bins=200

if True:
    hist, bin_edges = numpy.histogram(time_residuals_ppc[(dca_ppc>3.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(1.0,0.0,0.0), label=r"ppc ($\mathrm{pca} > 3 \,\mathrm{m}$)")
    hist, bin_edges = numpy.histogram(time_residuals_clsim[(dca_clsim>3.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(0.0,0.0,1.0), label=r"clsim ($\mathrm{pca} > 3 \,\mathrm{m}$)")

if False:
    hist, bin_edges = numpy.histogram(time_residuals_ppc[(dca_ppc<=5.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(1.0,0.0,0.0), label=r"ppc ($\mathrm{pca} \leq 5 \,\mathrm{m}$)")
    hist, bin_edges = numpy.histogram(time_residuals_clsim[(dca_clsim<=5.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(0.0,0.0,1.0), label=r"clsim ($\mathrm{pca} \leq 5 \,\mathrm{m}$)")

    hist, bin_edges = numpy.histogram(time_residuals_ppc[(dca_ppc>5.) & (dca_ppc<=10.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(1.0,0.2,0.2), label=r"ppc ($5\,\mathrm{m} < \mathrm{pca} \leq 10\,\mathrm{m}$)")
    hist, bin_edges = numpy.histogram(time_residuals_clsim[(dca_clsim>5.) & (dca_clsim<=10.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(0.2,0.2,1.0), label=r"clsim ($5\,\mathrm{m} < \mathrm{pca} \leq 10\,\mathrm{m}$)")

    hist, bin_edges = numpy.histogram(time_residuals_ppc[(dca_ppc>10.) & (dca_ppc<=20.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(1.0,0.4,0.4), label=r"ppc ($10\,\mathrm{m} < \mathrm{pca} \leq 20\,\mathrm{m}$)")
    hist, bin_edges = numpy.histogram(time_residuals_clsim[(dca_clsim>10.) & (dca_clsim<=20.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(0.4,0.4,1.0), label=r"clsim ($10\,\mathrm{m} < \mathrm{pca} \leq 20\,\mathrm{m}$)")

    hist, bin_edges = numpy.histogram(time_residuals_ppc[(dca_ppc>20.) & (dca_ppc<=50.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(1.0,0.6,0.6), label=r"ppc ($20\,\mathrm{m} < \mathrm{pca} \leq 50\,\mathrm{m}$)")
    hist, bin_edges = numpy.histogram(time_residuals_clsim[(dca_clsim>20.) & (dca_clsim<=50.)], bins=num_bins, range=the_range)
    bx.semilogy(bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2., hist, color=(0.6,0.6,1.0), label=r"clsim ($20\,\mathrm{m} < \mathrm{pca} \leq 50\,\mathrm{m}$)")




bx.legend(loc='upper right')
bx.grid(True)
bx.set_xlabel(r"$\delta t$")
bx.set_ylabel("number of hits")
bx.set_xlim(the_range[0],the_range[1])



print("saving as {0}".format(outfile))
pylab.savefig(outfile, transparent=False)



