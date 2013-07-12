#!/usr/bin/env python

from __future__ import print_function

from icecube import icetray, dataclasses, phys_services, clsim
from I3Tray import I3Units

import numpy
import scipy
import scipy.interpolate
import scipy.integrate

def setup_converter(useGeant4=False):
    # make a converter
    if useGeant4:
        ppcConverter = clsim.I3CLSimLightSourceToStepConverterGeant4()
    else:
        ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)

    # initialize it
    randomGen = phys_services.I3SPRNGRandomService(
        seed = 123456,
        nstreams = 10000,
        streamnum = 1)
    mediumProperties = clsim.MakeIceCubeMediumProperties()

    #DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    #RadiusOverSizeFactor = 5.
    #domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*RadiusOverSizeFactor)
    domAcceptance = clsim.I3CLSimFunctionConstant(1.)

    # lets set it up
    ppcConverter.SetMediumProperties(mediumProperties)
    ppcConverter.SetRandomService(randomGen)
    ppcConverter.SetWlenBias(domAcceptance)

    ppcConverter.SetMaxBunchSize(10240)
    ppcConverter.SetBunchSizeGranularity(1)

    ppcConverter.Initialize()

    return ppcConverter

def gen_steps(particle, converter, copies=1):
    
    # insert the requested number of copies of the particle
    for i in range(copies):
        # make a light source from the particle
        lightSource = clsim.I3CLSimLightSource(particle)
        
        # put it in the queue
        converter.EnqueueLightSource(lightSource, i)

    # tell the converter that we want all the results now
    converter.EnqueueBarrier()
    
    # retrieve all results
    steps = []
    while True:
        #barrierReset=False
        #stepSeries = converter.GetConversionResultWithBarrierInfo(barrierReset)
        
        stepSeries = converter.GetConversionResult()
        barrierReset = not converter.BarrierActive()
        
        for step in stepSeries:
            steps.append(step)
            #print step
        
        # get out of the loop if the barrier has been reset
        if barrierReset:
            break 
    
    return steps

# set up the Geant4 environment
clsim.AutoSetGeant4Environment()

# set up converter
ppcConverter = setup_converter()

def generate_stuff_at_energy(energy, iterations=100, copies=10):
    # make a particle
    p = dataclasses.I3Particle()
    p.pos = dataclasses.I3Position(0.,0.,-100.*I3Units.m)
    p.dir = dataclasses.I3Direction(1.,0.,0.)
    p.time = 0.
    p.energy = energy
    p.shape = dataclasses.I3Particle.ParticleShape.Cascade
    p.type = dataclasses.I3Particle.ParticleType.EMinus
    p.length = 0.
    p.location_type = dataclasses.I3Particle.LocationType.InIce

    weights = []
    xPos = []
    dirCosines = []
    dirAngles = []
    totalNumPhotons=0

    gran=iterations/10
    if gran <= 0: gran=1
    for it in range(iterations):
        if it%gran==0: print(it)
    
        # generate steps
        steps = gen_steps(p, ppcConverter, copies=copies)

        # calculate sums and other things
        for step in steps:
            totalNumPhotons += step.num
            weights.append(float(step.num))
            xPos.append(step.x)
            
            cosValue = step.dir.x*p.dir.x + step.dir.y*p.dir.y + step.dir.z*p.dir.z
            if cosValue < -1.: cosValue=-1.
            if cosValue > 1.: cosValue=1.
            dirCosines.append(cosValue)
            
            
    dirAngles = numpy.arccos(dirCosines) * 180./numpy.pi
    weights = numpy.array(weights)/float(iterations*copies)/(energy/I3Units.GeV)

    # some output
    print("generated", len(steps), "steps")
    print("with a total of", float(totalNumPhotons)/float(copies*iterations), " photons per cascade")


    x_pos_hist_data, x_pos_hist_edges = numpy.histogram(xPos, weights=weights, range=(-1.,8.), bins=500)
    x_pos_hist_centers = (x_pos_hist_edges[1:]+x_pos_hist_edges[:-1])/2.
    x_pos_hist_data = x_pos_hist_data / float(totalNumPhotons)


    dir_cos_hist_data, dir_cos_hist_edges = numpy.histogram(dirCosines, weights=weights, range=(-1.,1.), bins=200)
    dir_cos_hist_centers = (dir_cos_hist_edges[1:]+dir_cos_hist_edges[:-1])/2.
    dir_cos_hist_data = dir_cos_hist_data / float(totalNumPhotons)


    dir_ang_hist_data, dir_ang_hist_edges = numpy.histogram(dirAngles, weights=weights, range=(0.,180.), bins=200)
    dir_ang_hist_centers = (dir_ang_hist_edges[1:]+dir_ang_hist_edges[:-1])/2.
    dir_ang_hist_data = dir_ang_hist_data / float(totalNumPhotons)


    result = dict()
    result["total_num_photons"] = float(totalNumPhotons)/float(copies*iterations)
    result["total_num_steps"] = steps
    result["x_pos_hist_centers"] = x_pos_hist_centers
    result["x_pos_hist_data"] = x_pos_hist_data
    result["dir_cos_hist_centers"] = dir_cos_hist_centers
    result["dir_cos_hist_data"] = dir_cos_hist_data
    result["dir_ang_hist_centers"] = dir_ang_hist_centers
    result["dir_ang_hist_data"] = dir_ang_hist_data
    
    return result


energies =   [0.1*I3Units.GeV, 1.*I3Units.GeV, 10.*I3Units.GeV, 100.*I3Units.GeV, 1.*I3Units.TeV]
iterations = [100,             100,            10,              10,               1]
copies =     [500,             50,             50,              5,                5]

#energies = [0.1*I3Units.GeV, 1.*I3Units.GeV]

results = dict()
for i, energy in enumerate(energies):
    print("current energy: %fGeV" % (energy/I3Units.GeV))
    results[energy] = generate_stuff_at_energy(energy, iterations=iterations[i], copies=copies[i])

###############################


import matplotlib
matplotlib.use("PDF")

fig_size = [8.3,11.7] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': False,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
#matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

import pylab

def addAnnotationToPlot(plot, text, loc=1, size=8., rotation=0.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size, rotation=rotation),
                      frameon=True,
                      loc=loc, # 1=='upper right'
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)



fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)

for energy, result in sorted(results.items()):
    ax.plot(result["x_pos_hist_centers"], result["x_pos_hist_data"], label=r"E=%fGeV" % (energy/I3Units.GeV))

for energy, result in sorted(results.items()):
    bx.semilogy(result["dir_ang_hist_centers"], result["dir_ang_hist_data"], label=r"E=%fGeV" % (energy/I3Units.GeV))

energies=[]
num_photons=[]
for energy, result in sorted(results.items()):
    energies.append(energy)
    num_photons.append(float(result["total_num_photons"]))
energies=numpy.array(energies)
num_photons=numpy.array(num_photons)
cx.scatter(energies/I3Units.GeV, num_photons)
cx.loglog(energies/I3Units.GeV, num_photons)


ax.set_ylim(0.,0.000008)
ax.legend()
ax.grid(True)
#ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
#ax.set_ylabel("DOM acceptance")

bx.set_ylim(1e-10,1e-4)
bx.legend()
bx.grid(True)

cx.set_ylim(1e4,1e9)
cx.legend()
cx.grid(True)


pylab.savefig("cascadeGeneratorPlots.pdf")




