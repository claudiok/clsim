from __future__ import print_function
from I3Tray import I3Units
from icecube import icetray, dataclasses, dataio
from icecube import clsim

import math
import numpy

#import pylab
#import scipy
#import scipy.interpolate
#import scipy.integrate

import matplotlib
matplotlib.use("PDF")
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

def addAnnotationToPlot(plot, text, loc=1, size=8., rotation=0.):
    at = AnchoredText(text,
                      prop=dict(size=size, rotation=rotation),
                      frameon=True,
                      loc=loc, # 1=='upper right'
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)

def plotDOMs(ax, fig, frame, radius, minZ, maxZ, pi0Only=False, eplusOnly=False):
    #propagatedPhotons = frame["PropagatedPhotons"]
    #isHits = False
    
    propagatedPhotons = frame["MCHitSeriesMap"]
    isHits = True
    
    geometry = frame["I3Geometry"]
    mcTree = frame["I3MCTree"]
    
    primaries = mcTree.GetPrimaries()
    if len(primaries) > 1: raise RuntimeError("more than one primary..")
    
    if primaries[0].GetType() != dataclasses.I3Particle.ParticleType.PPlus:
        raise RuntimeError("Not a proton primary!")
    daughters = mcTree.GetDaughters(primaries[0])
    if len(daughters) != 2: raise RuntimeError("NumDaughters != 2")
    
    pi0Particle=None
    eplusParticle=None
    for daughter in daughters:
        if daughter.GetType() == dataclasses.I3Particle.ParticleType.Pi0:
            pi0Particle = daughter
        elif daughter.GetType() == dataclasses.I3Particle.ParticleType.EPlus or daughter.GetType() == dataclasses.I3Particle.ParticleType.MuPlus:
            eplusParticle = daughter
    if pi0Particle is None: raise RuntimeError("No pi0 found!")
    if eplusParticle is None: raise RuntimeError("No e+ found!")
    
    
    
    if False:
        firstTime=None
        lastTime=None
        
        for key, photonList in propagatedPhotons:
            for photon in photonList:
                if photon.numScattered > 0: continue
                if pi0Only and (photon.particleMajorID!=pi0Particle.major_id or photon.particleMinorID!=pi0Particle.minor_id): continue
                if eplusOnly and (photon.particleMajorID!=eplusParticle.major_id or photon.particleMinorID!=eplusParticle.minor_id): continue

                if firstTime is None:
                    firstTime=photon.time
                    lastTime=photon.time
                    continue
                if photon.time < firstTime: firstTime=photon.time
                if photon.time > lastTime: lastTime=photon.time

        if (firstTime is None) or (firstTime==lastTime):
            for key, photonList in propagatedPhotons:
                for photon in photonList:
                    if pi0Only and (photon.particleMajorID!=pi0Particle.major_id or photon.particleMinorID!=pi0Particle.minor_id): continue
                    if eplusOnly and (photon.particleMajorID!=eplusParticle.major_id or photon.particleMinorID!=eplusParticle.minor_id): continue

                    if firstTime is None:
                        firstTime=photon.time
                        lastTime=photon.time
                        continue
                    if photon.time < firstTime: firstTime=photon.time
                    if photon.time > lastTime: lastTime=photon.time

        
    
    plotPositionsX = []
    plotPositionsY = []
    plotPositionsColor = []
    plotPositionsSize = []

    numPhotons=0
    numDOMs=0
    
    #print "firstTime =", firstTime
    #    print " lastTime =", lastTime
    
    for key, photonList in propagatedPhotons:
        photonsInDOM=0
        currentPhoton=None
        for photon in photonList:
            #if photon.numScattered > 1: continue
            if isHits:
                if pi0Only and (photon.particleMajorID!=pi0Particle.major_id or photon.particleMinorID!=pi0Particle.minor_id): continue
                if eplusOnly and (photon.particleMajorID!=eplusParticle.major_id or photon.particleMinorID!=eplusParticle.minor_id): continue
            else:
                if pi0Only and (photon.particle_major_id!=pi0Particle.major_id or photon.particle_minor_id!=pi0Particle.minor_id): continue
                if eplusOnly and (photon.particle_major_id!=eplusParticle.major_id or photon.particle_minor_id!=eplusParticle.minor_id): continue

            if currentPhoton is None:
                currentPhoton=photon
            elif photon.time < currentPhoton.time:
                currentPhoton=photon
            numPhotons+=1
            photonsInDOM+=1
        if currentPhoton is None: continue
        
        numDOMs+=1
        
        posInTimeInterval = (currentPhoton.time-firstTime)/(lastTime-firstTime)
        
        #if posInTimeInterval < 0.: posInTimeInterval=0.
        if posInTimeInterval > 1.: posInTimeInterval=1.
        #if posInTimeInterval < firstTime: posInTimeInterval=firstTime
        #if posInTimeInterval > lastTime: posInTimeInterval=lastTime
                        
        omgeo = geometry.omgeo[key]
        
        zCylPos = omgeo.position.Z
        phiCylPos = math.atan2(omgeo.position.Y,omgeo.position.X)
        
        plotPositionsX.append(phiCylPos)
        plotPositionsY.append(zCylPos)
        #plotPositionsColor.append(posInTimeInterval)
        plotPositionsColor.append(currentPhoton.time/I3Units.ns)
        
        dotSize = float(photonsInDOM)*3.
        if dotSize > 20.: dotSize=20.
        plotPositionsSize.append(dotSize)
    
    plotPositionsX = numpy.array(plotPositionsX)
    plotPositionsY = numpy.array(plotPositionsY)
    plotPositionsSize = numpy.array(plotPositionsSize)
    
    l = ax.scatter((plotPositionsX+math.pi)*radius/I3Units.m, plotPositionsY/I3Units.m, c=plotPositionsColor, s=plotPositionsSize, edgecolors="none", norm=matplotlib.colors.Normalize(0.,250.))
    
    cbar = fig.colorbar(l, ax=ax)
    

    statistics = frame["CLSimStatistics"]
    if statistics is not None:
        massPi0 = 134.9766*I3Units.MeV
        massPositron = 0.510998910*I3Units.MeV
        
        ETotPi0 = pi0Particle.GetEnergy()+massPi0
        ETotPositron = eplusParticle.GetEnergy()+massPositron
        
        numPhotonsAtDOMs = statistics.GetTotalSumOfWeightsPhotonsAtDOMs()
        numPhotonsGen = statistics.GetTotalSumOfWeightsPhotonsGenerated()

        numPhotonsAtDOMs_pi0 = statistics.GetSumOfWeightsPhotonsAtDOMsForParticle(pi0Particle)
        numPhotonsGen_pi0 = statistics.GetSumOfWeightsPhotonsGeneratedForParticle(pi0Particle)

        numPhotonsAtDOMs_eplus = statistics.GetSumOfWeightsPhotonsAtDOMsForParticle(eplusParticle)
        numPhotonsGen_eplus = statistics.GetSumOfWeightsPhotonsGeneratedForParticle(eplusParticle)

        addAnnotationToPlot(ax,
        r"\begin{flushright}" +
        r"total: photons generated: %g, photons @ DOMs: %g, ratio: %.1f\%%\\" % (numPhotonsGen, numPhotonsAtDOMs, 100.*numPhotonsAtDOMs/numPhotonsGen) +
        r"$\pi^0$ ($E_\mathrm{tot}=%.2f\,\mathrm{MeV}$): photons generated ($N_\mathrm{gen}$): %g, photons @ DOMs: %g, ratio: %.1f\%%\\" % (ETotPi0/I3Units.MeV, numPhotonsGen_pi0, numPhotonsAtDOMs_pi0, 100.*numPhotonsAtDOMs_pi0/numPhotonsGen_pi0) +
        r"$\mathrm{e}^+$ ($E_\mathrm{tot}=%.2f\,\mathrm{MeV}$): photons generated ($N_\mathrm{gen}$): %g, photons @ DOMs: %g, ratio: %.1f\%%" % (ETotPositron/I3Units.MeV, numPhotonsGen_eplus, numPhotonsAtDOMs_eplus, 100.*numPhotonsAtDOMs_eplus/numPhotonsGen_eplus) +
        r"\end{flushright}")
    
    
        photonsGeneratedPerEnergy_pi0 = float(numPhotonsGen_pi0)/ETotPi0
        photonsGeneratedPerEnergy_eplus = float(numPhotonsGen_eplus)/ETotPositron

        addAnnotationToPlot(ax,
        r"\begin{flushleft}" +
        r"$\pi^0$: $N_\mathrm{gen} / E_\mathrm{tot} = %g\,\mathrm{MeV}^{-1}$\\" % (photonsGeneratedPerEnergy_pi0/(1./I3Units.MeV)) +
        r"$\mathrm{e}^+$: $N_\mathrm{gen} / E_\mathrm{tot} = %g\,\mathrm{MeV}^{-1}$\\" % (photonsGeneratedPerEnergy_eplus/(1./I3Units.MeV)) +
        r"\end{flushleft}", loc=2)
    
    
    ax.fill_between([-1.,2.*math.pi*(radius/I3Units.m)+1.], [maxZ, maxZ], [maxZ+200., maxZ+200.], linewidth=2., color='k', alpha=0.5)
    
    ax.set_xlim(0.,2.*math.pi*radius/I3Units.m)
    ax.set_ylim(minZ,maxZ+(maxZ-minZ)*0.1) # some extra space at the top
    #ax.legend()
    ax.grid(True)
    ax.set_ylabel(r"$z$ [$\mathrm{m}$]")
    ax.set_xlabel(r"$R \times \varphi$ [$\mathrm{m}$] ($R=%g\,\mathrm{m}$)" % (radius/I3Units.m))
    ax.set_aspect('equal', 'box')
    #ax.set_aspect('equal', 'datalim')
    cbar.set_label(r"absolute time $t_\mathrm{abs}$ [$\mathrm{ns}$]")
    
    
    
    print("numPhotons={0}".format(numPhotons))
    print("   numDOMs={0}".format(numDOMs))


class eventPlotter(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("filename", "output pdf filename", "")
        self.AddParameter("pi0Only", "only plot hits from pi0", False)
        self.AddParameter("eplusOnly", "only plot hits from e+", False)
        self.AddOutBox("OutBox")
        
        self.plotNum=0
        
    def Configure(self):
        self.filename = self.GetParameter("filename")
        self.pi0Only = self.GetParameter("pi0Only")
        self.eplusOnly = self.GetParameter("eplusOnly")
        
        from matplotlib.backends.backend_pdf import PdfPages
        self.pdffile = PdfPages(self.filename)
        
        fig_size = [11.7,8.3] # din A4
        params = {'backend': 'pdf',
                'axes.labelsize': 12,
                'text.fontsize': 12,
                'legend.fontsize': 12,
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'text.usetex': True,
                'figure.figsize': fig_size}
        matplotlib.rcParams.update(params)
        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
        
        
    def Physics(self, frame):
        eventHeader = frame["I3EventHeader"]
        
        maxR=0.
        meanR=0.
        entries=0
        maxZ=None
        minZ=None
        geometry = frame["I3Geometry"]
        stringnumset = set()
        domZPosPrev = None
        domZSpacing = None
        numDOMs=0

        for key, omgeo in geometry.omgeo:
            if key.GetString() <= 86: continue
            
            stringnumset.add(key.GetString())
            numDOMs+=1

            if domZPosPrev is None:
                domZPosPrev = omgeo.position.Z
            elif domZSpacing is None:
                domZSpacing = omgeo.position.Z-domZPosPrev

            if maxZ is None:
                maxZ=omgeo.position.Z
                minZ=omgeo.position.Z
            else:
                if omgeo.position.Z > maxZ: maxZ=omgeo.position.Z
                if omgeo.position.Z < minZ: minZ=omgeo.position.Z
               
            radius = math.sqrt(omgeo.position.X**2. + omgeo.position.Y**2.)
            if radius > maxR: maxR=radius
            meanR += radius
            entries += 1
        meanR /= float(entries)

        print("    radius={0}m".format(maxR/I3Units.m))

        #minZ=-500.
        #maxZ=-100.
        #meanR=33.
        #maxR=33.

        fig = matplotlib.pyplot.figure(facecolor='white')
        fig.subplots_adjust(left=0.09, bottom=0.055, top=0.93, right=0.98)

        ax = fig.add_subplot(1, 1, 1)

        plotDOMs(ax, fig, frame, maxR, minZ, maxZ, pi0Only=self.pi0Only, eplusOnly=self.eplusOnly)

        extraString=""
        if self.pi0Only: extraString += r"$\pi^0$ only; "
        if self.eplusOnly: extraString += r"$\mathrm{e}^-$ only; "

        fig.suptitle(r"\begin{center}" + extraString +
                     r"(run %u event %u) " % (eventHeader.RunID, eventHeader.EventID) +
                     r"perfect photon counting (all photons $\lambda \in \left[ 265\,\mathrm{nm} ; 675\,\mathrm{nm} \right]$)\\" +
                     r"IceCube coordinates, ref. depth ($z=0$) is $%g\,\mathrm{m}$;"%(1948.07) +
                     r"$N_\mathrm{string}=%u$; $d_{\mathrm{DOM};z}=%g\,\mathrm{m}$; $N_\mathrm{DOM}=%u$ " % (len(stringnumset), domZSpacing/I3Units.m, numDOMs) +
                     r"\end{center}", fontsize=12)

        #matplotlib.pyplot.savefig(self.filename, transparent=False)
        print("saving plot (page %u)." % (self.plotNum))
        fig.savefig(self.pdffile, format='pdf')
        self.plotNum+=1
        
        
        self.PushFrame(frame)
        
        
    def Finish(self):
        print("closing pdf file...")
        self.pdffile.close()
        print("finished!")
