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

import matplotlib.projections
import matplotlib.cm


def getPhaseRefIndex(wavelength):
    x = wavelength/1000.# wavelength in micrometer
    return 1.55749 - 1.57988*x + 3.99993*x**2. - 4.68271*x**3. + 2.09354*x**4.
def getGroupRefIndex(wavelength):
    np = getPhaseRefIndex(wavelength)
    x = wavelength/1000.# wavelength in micrometer
    return np * (1. + 0.227106 - 0.954648*x + 1.42568*x**2. - 0.711832*x**3.)


def addAnnotationToPlot(plot, text, loc=1, size=8., rotation=0., bbox_transform=None, bbox_to_anchor=None):
    at = AnchoredText(text,
                      prop=dict(size=size, rotation=rotation),
                      frameon=True,
                      loc=loc, # 1=='upper right',
                      bbox_transform=bbox_transform,
                      bbox_to_anchor=bbox_to_anchor
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)

def plotDOMs(fig, frame, meanPosX, meanPosY, radius, minZ, maxZ, pi0Only=False, eplusOnly=False, subdetectorName="MICA"):
    #propagatedPhotons = frame["PhotonsOnDOMs"]
    #isHits = False
    
    propagatedPhotons = frame["MCHitSeriesMap_MICA_mDOM"]
    isHits = True
    
    moduleGeoMap = frame["I3ModuleGeoMap"]
    subdetectors = frame["Subdetectors"]
    
    mcTree = frame["I3MCTree"]
    
    primaries = mcTree.primaries
    if len(primaries) > 1: raise RuntimeError("more than one primary..")
    
    if primaries[0].type != dataclasses.I3Particle.ParticleType.PPlus:
        raise RuntimeError("Not a proton primary!")
    daughters = mcTree.get_daughters(primaries[0])
    if len(daughters) != 2: raise RuntimeError("NumDaughters != 2")
    
    pi0Particle=None
    eplusParticle=None
    for daughter in daughters:
        if daughter.type == dataclasses.I3Particle.ParticleType.Pi0:
            pi0Particle = daughter
        elif daughter.type == dataclasses.I3Particle.ParticleType.EPlus or daughter.type == dataclasses.I3Particle.ParticleType.MuPlus:
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
    
    protonR = math.sqrt((primaries[0].pos.x-meanPosX)**2 + (primaries[0].pos.y-meanPosY)**2)
    protonZ = primaries[0].pos.z
    
    #print "firstTime =", firstTime
    #    print " lastTime =", lastTime
    
    for key, photonList in propagatedPhotons:
        if isHits:
            if subdetectors[dataclasses.ModuleKey(key.string, key.om)] != subdetectorName: continue
        else:
            if subdetectors[key] != subdetectorName: continue
        
        photonsInDOM=0
        currentPhoton=None
        for photon in photonList:
            #if photon.numScattered > 1: continue
            if not isHits:
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
        
        #posInTimeInterval = (currentPhoton.time-firstTime)/(lastTime-firstTime)
        
        #if posInTimeInterval < 0.: posInTimeInterval=0.
        #if posInTimeInterval > 1.: posInTimeInterval=1.
        #if posInTimeInterval < firstTime: posInTimeInterval=firstTime
        #if posInTimeInterval > lastTime: posInTimeInterval=lastTime
                        
        modulegeo = moduleGeoMap[dataclasses.ModuleKey(key.string, key.om)]
        
        relativePosX = modulegeo.pos.x - primaries[0].pos.x
        relativePosY = modulegeo.pos.y - primaries[0].pos.y
        relativePosZ = modulegeo.pos.z - primaries[0].pos.z        

        r = math.sqrt(relativePosX**2 + relativePosY**2 + relativePosZ**2)
        phi = math.atan2(relativePosY,relativePosX)
        theta = math.acos(relativePosZ/r)-math.pi/2.

        plotTime = currentPhoton.time/I3Units.ns
        
        if isHits:
            expectedTime = r/(dataclasses.I3Constants.c/getGroupRefIndex(400.))
        else:
            expectedTime = r/(dataclasses.I3Constants.c/getGroupRefIndex(currentPhoton.wavelength/I3Units.nanometer))
        plotTime = (currentPhoton.time-expectedTime)/I3Units.ns
        #print plotTime
        #if plotTime > 10.: continue
        
        plotPositionsX.append(phi)
        plotPositionsY.append(theta)
        plotPositionsColor.append(plotTime)
        
        dotSize = float(photonsInDOM)*3.
        if dotSize > 200.: dotSize=200.
        plotPositionsSize.append(dotSize)
    
    plotPositionsX = numpy.array(plotPositionsX)
    plotPositionsY = numpy.array(plotPositionsY)
    plotPositionsSize = numpy.array(plotPositionsSize)
    
    
    fig.subplots_adjust(left=0.02, bottom=0.055, top=0.93, right=0.98)
    
    
    mc_pi0_phi = math.atan2((pi0Particle.dir.y-meanPosY), (pi0Particle.dir.x-meanPosX))
    mc_pi0_theta = math.acos(pi0Particle.dir.z)-math.pi/2.
    mc_eplus_phi = math.atan2((eplusParticle.dir.y-meanPosY), (eplusParticle.dir.x-meanPosX))
    mc_eplus_theta = math.acos(eplusParticle.dir.z)-math.pi/2.

    ax = fig.add_subplot(1, 2, 1, projection='lambert', center_longitude=mc_pi0_phi, center_latitude=mc_pi0_theta)
    bx = fig.add_subplot(1, 2, 2, projection='lambert', center_longitude=mc_eplus_phi, center_latitude=mc_eplus_theta)
    
    startColor=(0.1,0.1,0.1)
    
    midPointStart=0.3
    midStartColor=(1.0,0.0,0.0)
    midEndColor=(0.4,0.0,0.0)
    midPointEnd=0.7
    
    endColor=(0.8,0.8,0.8)
    
    #colorMap = matplotlib.cm.Set1
    cdict = {'red': ((0.0, 0.0, startColor[0]),
                     (midPointStart, startColor[0] + (endColor[0]-startColor[0])*midPointStart, midStartColor[0]),
                     (midPointEnd,   midEndColor[0], startColor[0] + (endColor[0]-startColor[0])*midPointEnd),
                     (1.0, endColor[0], 1.0)),
           'green': ((0.0, 0.0, startColor[1]),
                     (midPointStart, startColor[1] + (endColor[1]-startColor[1])*midPointStart, midStartColor[1]),
                     (midPointEnd,   midEndColor[1], startColor[1] + (endColor[1]-startColor[1])*midPointEnd),
                     (1.0, endColor[1], 1.0)),
            'blue': ((0.0, 0.0, startColor[2]),
                     (midPointStart, startColor[2] + (endColor[2]-startColor[2])*midPointStart, midStartColor[2]),
                     (midPointEnd,   midEndColor[2], startColor[2] + (endColor[2]-startColor[2])*midPointEnd),
                     (1.0, endColor[2], 1.0))}
    colorMap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #colorMap = None
    
    l = ax.scatter(plotPositionsX, plotPositionsY, c=plotPositionsColor, s=plotPositionsSize, cmap=colorMap, edgecolors="none", norm=matplotlib.colors.Normalize(-10.,35.))
    cbar_ax = fig.add_axes([0.02, 0.1, 0.98-0.02, 0.05])
    cbar = fig.colorbar(l, cax=cbar_ax, orientation='horizontal')
    
    l = bx.scatter(plotPositionsX, plotPositionsY, c=plotPositionsColor, s=plotPositionsSize, cmap=colorMap, edgecolors="none", norm=matplotlib.colors.Normalize(-10.,35.))

    #label_e  = ax.scatter([mc_eplus_phi, mc_eplus_phi+1.], [mc_eplus_theta, mc_eplus_theta+0.0001], s=[50., 0.], color='b', label=r"$\mathrm{e}^{+}$")
    label_pi = ax.scatter([mc_pi0_phi, mc_pi0_phi+1.],     [mc_pi0_theta, mc_pi0_theta+0.0001],     s=[50., 0.], color='y', label=r"$\pi^{0}$")
    #label_pi = bx.scatter([mc_pi0_phi, mc_pi0_phi+1.],     [mc_pi0_theta, mc_pi0_theta+0.0001],     s=[50., 0.], color='y', label=r"$\pi^{0}$")
    label_e  = bx.scatter([mc_eplus_phi, mc_eplus_phi+1.], [mc_eplus_theta, mc_eplus_theta+0.0001], s=[50., 0.], color='b', label=r"$\mathrm{e}^{+}$")

    #print mc_pi0_theta
    #label_pi = bx.scatter([mc_pi0_phi], [1.], s=[50.], color='b', label=r"$\mathrm{e}^{+}$")
    #label_e  = bx.scatter([mc_eplus_phi], [mc_eplus_theta], s=[50.], color='b', label=r"$\mathrm{e}^{+}$")

        
    addAnnotationToPlot(ax,
        r"$r_\mathrm{proton}=%f\;\mathrm{m} \;;\; z_\mathrm{proton}=%f\;\mathrm{m}$" % (protonR/I3Units.m, protonZ/I3Units.m),
        loc=9, bbox_transform=fig.transFigure, bbox_to_anchor=(0., 0., 1., .85), size=10.)
    
    statistics = frame["CLSimStatistics"]
    if statistics is not None:
        massPi0 = 134.9766*I3Units.MeV
        massPositron = 0.510998910*I3Units.MeV
    
        ETotPi0 = pi0Particle.energy+massPi0
        ETotPositron = eplusParticle.energy+massPositron
    
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
        r"\end{flushright}",
        bbox_transform=fig.transFigure, bbox_to_anchor=(0.01, 0., 0.98, .93),
        loc=1)


        photonsGeneratedPerEnergy_pi0 = float(numPhotonsGen_pi0)/ETotPi0
        photonsGeneratedPerEnergy_eplus = float(numPhotonsGen_eplus)/ETotPositron

        addAnnotationToPlot(ax,
        r"\begin{flushleft}" +
        r"$\pi^0$: $N_\mathrm{gen} / E_\mathrm{tot} = %g\,\mathrm{MeV}^{-1}$\\" % (photonsGeneratedPerEnergy_pi0/(1./I3Units.MeV)) +
        r"$\mathrm{e}^+$: $N_\mathrm{gen} / E_\mathrm{tot} = %g\,\mathrm{MeV}^{-1}$\\" % (photonsGeneratedPerEnergy_eplus/(1./I3Units.MeV)) +
        r"\end{flushleft}",
        bbox_transform=fig.transFigure, bbox_to_anchor=(0.01, 0., 0.98, .93),
        loc=2)
    

    ax.legend([label_pi], [r"$\pi^0$"], scatterpoints=1, loc='upper right')
    bx.legend([label_e], [r"$\mathrm{e}^{+}$"], scatterpoints=1, loc='upper left')
    
    ax.grid(True)
    bx.grid(True)

    cbar.set_label(r"time delay $t_\mathrm{detected} - t_\mathrm{expected}$ [$\mathrm{ns}$]")
    
    
    
    print("numPhotons={0}".format(numPhotons))
    print("   numDOMs={0}".format(numDOMs))

class eventPlotter2(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("filename", "output pdf filename", "")
        self.AddParameter("pi0Only", "only plot hits from pi0", False)
        self.AddParameter("eplusOnly", "only plot hits from e+", False)
        self.AddParameter("SubdetectorName", "Only plot hits in this subdetector", "MICA")
        self.AddOutBox("OutBox")
        
        self.plotNum=0
        
    def Configure(self):
        self.filename = self.GetParameter("filename")
        self.pi0Only = self.GetParameter("pi0Only")
        self.eplusOnly = self.GetParameter("eplusOnly")
        self.subdetectorName = self.GetParameter("SubdetectorName")
        
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
        
        
    def DAQ(self, frame):
        eventHeader = frame["I3EventHeader"]
        
        maxR=0.
        meanR=0.
        entries=0
        maxZ=None
        minZ=None
        moduleGeoMap = frame["I3ModuleGeoMap"]
        subdetectors = frame["Subdetectors"]
        stringnumset = set()
        domZPosPrev = None
        domZSpacing = None
        numDOMs=0

        meanPosX = 0.
        meanPosY = 0.
        counter = 0
        for key, modulegeo in moduleGeoMap:
            if subdetectors[key] != self.subdetectorName: continue

            meanPosX += modulegeo.pos.x
            meanPosY += modulegeo.pos.y
            counter += 1
        meanPosX /= float(counter)
        meanPosY /= float(counter)

        for key, modulegeo in moduleGeoMap:
            if subdetectors[key] != self.subdetectorName: continue
            
            stringnumset.add(key.string)
            numDOMs+=1

            if domZPosPrev is None:
                domZPosPrev = modulegeo.pos.z
            elif domZSpacing is None:
                domZSpacing = modulegeo.pos.z-domZPosPrev

            if maxZ is None:
                maxZ=modulegeo.pos.z
                minZ=modulegeo.pos.z
            else:
                if modulegeo.pos.z > maxZ: maxZ=modulegeo.pos.z
                if modulegeo.pos.z < minZ: minZ=modulegeo.pos.z
               
            radius = math.sqrt((modulegeo.pos.x-meanPosX)**2. + (modulegeo.pos.y-meanPosY)**2.)
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

        #ax = fig.add_subplot(1, 1, 1, projection='lambert', center_longitude=math.pi/4., center_latitude=math.pi/4.)
        #ax = fig.add_subplot(1, 1, 1, projection='lambert')

        plotDOMs(fig, frame, meanPosX, meanPosY, maxR, minZ, maxZ, self.subdetectorName)

        extraString=""
        if self.pi0Only: extraString += r"$\pi^0$ only; "
        if self.eplusOnly: extraString += r"$\mathrm{e}^-$ only; "

        fig.suptitle(r"\begin{center}" + extraString +
                     r"(run %u event %u) " % (eventHeader.run_id, eventHeader.event_id) +
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
