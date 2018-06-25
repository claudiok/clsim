#!/usr/bin/env python
from os.path import expandvars
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="MICA_proton_decay.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-g", "--gcd",default="geometry_MICA.i3",
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=54321,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

###############################

from I3Tray import *

from icecube import icetray, dataclasses, dataio, phys_services, sim_services
from icecube import clsim

from proton_decay_generator import I3SimpleProtonDecayGenerator

import math, numpy

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray = I3Tray()

# make a stream of events
tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)
tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=options.RUNNUMBER,
               EventID=1,
               IncrementEventID=True)
tray.AddModule("Delete", "cleanup", Keys=["MCTimeIncEventID"])

# decay some protons
tray.AddModule(I3SimpleProtonDecayGenerator, "I3SimpleProtonDecayGenerator",
               #posZRange=(-327.*I3Units.m,-502.*I3Units.m),
               #radius=50.*I3Units.m,
               #centerX=0.*I3Units.m,
               #centerY=0.*I3Units.m,

               posZRange=(-500.*I3Units.m,-100.*I3Units.m),
               radius=60.*I3Units.m,
               centerX=45.771827*I3Units.m,
               centerY=-34.411053*I3Units.m,
   
               randomService=randomService,
               Streams=[icetray.I3Frame.DAQ])

# propagate particles, generate C'kov photons, track them and create hits in IceCube DOMs
tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
                ParallelEvents = 1,     # this assumes low energies, you will probably run out of memory if your energies are too high
                RandomService = randomService,
                UseGPUs=False,
                UseCPUs=True, 

                UseGeant4=True,         # use Geant4 (i.e. no parameterizations)
                DOMOversizeFactor=1.,
                UnweightedPhotons=True,
                PhotonSeriesName="PhotonsOnDOMs",

                MMCTrackListName=None,  # no MMC!
                IceModelLocation=expandvars("$I3_BUILD/ice-models/resources/models/spice_mie"),
                ExtraArgumentsToI3CLSimModule = dict(StatisticsName="CLSimStatistics")
                )


glassThickness = 14.*I3Units.mm
glassAbslen = [
                 0.00*I3Units.cm,  # 260nm
                 0.00*I3Units.cm,  # 270nm
                 0.00*I3Units.cm,  # 280nm
                 0.00*I3Units.cm,  # 290nm
                 
                 0.17*I3Units.cm,  # 300nm
                 0.39*I3Units.cm,  # 310nm
                 0.84*I3Units.cm,  # 320nm
                 1.82*I3Units.cm,  # 330nm
                 3.92*I3Units.cm,  # 340nm
                 8.41*I3Units.cm,  # 350nm
                18.09*I3Units.cm,  # 360nm
                27.21*I3Units.cm,  # 370nm
                19.23*I3Units.cm,  # 380nm
                61.84*I3Units.cm,  # 390nm
               128.04*I3Units.cm,  # 400nm
                81.25*I3Units.cm,  # 410nm
                73.02*I3Units.cm,  # 420nm
                77.30*I3Units.cm,  # 430nm
                65.66*I3Units.cm,  # 440nm
                81.63*I3Units.cm,  # 450nm
               109.23*I3Units.cm,  # 460nm
               116.08*I3Units.cm,  # 470nm
               113.90*I3Units.cm,  # 480nm
               118.86*I3Units.cm,  # 490nm
               126.55*I3Units.cm,  # 500nm
               139.70*I3Units.cm,  # 510nm
               145.68*I3Units.cm,  # 520nm
               150.88*I3Units.cm,  # 530nm
               151.80*I3Units.cm,  # 540nm
               147.16*I3Units.cm,  # 550nm
               142.40*I3Units.cm,  # 560nm
               138.27*I3Units.cm,  # 570nm
               134.58*I3Units.cm,  # 580nm
               135.64*I3Units.cm,  # 590nm
               142.87*I3Units.cm,  # 600nm
               148.37*I3Units.cm,  # 610nm
               
               148.37*I3Units.cm,  # 620nm
               148.37*I3Units.cm,  # 630nm
               148.37*I3Units.cm,  # 640nm
               148.37*I3Units.cm,  # 650nm
               148.37*I3Units.cm,  # 660nm
               148.37*I3Units.cm,  # 670nm
               148.37*I3Units.cm,  # 680nm
               148.37*I3Units.cm,  # 690nm
               
               ]
glassAbsLen = clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, glassAbslen)

gelAbsLen = [  0.00*I3Units.cm, # 260nm
               0.00*I3Units.cm, # 270nm
               0.00*I3Units.cm, # 280nm
               0.00*I3Units.cm, # 290nm
               
               0.00*I3Units.cm, # 300nm
               8.00*I3Units.cm, # 310nm
              15.60*I3Units.cm, # 320nm
              23.08*I3Units.cm, # 330nm
              30.49*I3Units.cm, # 340nm
              37.14*I3Units.cm, # 350nm
              41.88*I3Units.cm, # 360nm
              45.71*I3Units.cm, # 370nm
              48.96*I3Units.cm, # 380nm
              53.29*I3Units.cm, # 390nm
              56.64*I3Units.cm, # 400nm
              59.38*I3Units.cm, # 410nm
              62.53*I3Units.cm, # 420nm
              64.48*I3Units.cm, # 430nm
              66.91*I3Units.cm, # 440nm
              68.05*I3Units.cm, # 450nm
              72.31*I3Units.cm, # 460nm
              74.55*I3Units.cm, # 470nm
              76.48*I3Units.cm, # 480nm
              78.18*I3Units.cm, # 490nm
              81.08*I3Units.cm, # 500nm
              84.49*I3Units.cm, # 510nm
              85.88*I3Units.cm, # 520nm
              86.95*I3Units.cm, # 530nm
              90.10*I3Units.cm, # 540nm
              89.09*I3Units.cm, # 550nm
              94.36*I3Units.cm, # 560nm
              96.42*I3Units.cm, # 570nm
              96.90*I3Units.cm, # 580nm
              99.89*I3Units.cm, # 590nm
              99.94*I3Units.cm, # 600nm
             100.81*I3Units.cm, # 610nm

             100.81*I3Units.cm, # 620nm
             100.81*I3Units.cm, # 630nm
             100.81*I3Units.cm, # 640nm
             100.81*I3Units.cm, # 650nm
             100.81*I3Units.cm, # 660nm
             100.81*I3Units.cm, # 670nm
             100.81*I3Units.cm, # 680nm
             100.81*I3Units.cm, # 690nm
             ]
gelAbsLen = clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, gelAbsLen)

wavelengthAcceptance = [ 0.0 * 0.01, # 260nm
                         0.0 * 0.01, # 270nm
                         0.5 * 0.01, # 280nm
                         3.1 * 0.01, # 290nm
                         9.8 * 0.01, # 300nm
                        17.5 * 0.01, # 310nm
                        23.2 * 0.01, # 320nm
                        26.5 * 0.01, # 330nm
                        28.1 * 0.01, # 340nm
                        28.1 * 0.01, # 350nm
                        29.1 * 0.01, # 360nm
                        30.1 * 0.01, # 370nm
                        30.4 * 0.01, # 380nm
                        30.1 * 0.01, # 390nm
                        29.9 * 0.01, # 400nm
                        29.3 * 0.01, # 410nm
                        28.6 * 0.01, # 420nm
                        27.5 * 0.01, # 430nm
                        26.5 * 0.01, # 440nm
                        25.0 * 0.01, # 450nm
                        23.2 * 0.01, # 460nm
                        21.1 * 0.01, # 470nm
                        19.6 * 0.01, # 480nm
                        18.5 * 0.01, # 490nm
                        17.2 * 0.01, # 500nm
                        15.4 * 0.01, # 510nm
                        12.1 * 0.01, # 520nm
                         9.3 * 0.01, # 530nm
                         7.2 * 0.01, # 540nm
                         6.2 * 0.01, # 550nm
                         4.6 * 0.01, # 560nm
                         3.6 * 0.01, # 570nm
                         2.8 * 0.01, # 580nm
                         2.1 * 0.01, # 590nm
                         1.3 * 0.01, # 600nm
                         0.8 * 0.01, # 610nm
                         0.5 * 0.01, # 620nm
                         0.3 * 0.01, # 630nm
                         0.0 * 0.01, # 640nm
                         0.0 * 0.01, # 650nm
                         0.0 * 0.01, # 660nm
                         0.0 * 0.01, # 670nm
                         0.0 * 0.01, # 680nm
                         0.0 * 0.01, # 690nm
                         ]
wavelengthAcceptance = clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, wavelengthAcceptance)

cosines = [-1.,   -0.95, -0.9,  -0.85, -0.8,  -0.75, -0.7,  -0.65, -0.6 , -0.55, -0.5 , -0.45,
           -0.4,  -0.35, -0.3,  -0.25, -0.2,  -0.15, -0.1,  -0.05,  0.,    0.05]
acceptanceRatios = [ 1.56901889,   1.23410288,   1.17370277,   1.13023417,   1.1384318,
                     1.12859394,   1.09788943,   1.10675797,   1.11971921,   1.11214495,
                     1.15454716,   1.15321001,   1.2245022,    1.27309651,   1.30277673,
                     1.40193728,   1.52872749,   1.92281438,   1.94452437,   1.91083794,
                     1.92,         1.92] # regularize the acceptance ratios here..
                     #2.91631343,  29.62337522]
winstonRatioFunc = clsim.util.interpolate.interp1d(cosines,acceptanceRatios)

cos_bins = numpy.linspace(-1.,1.,1001)
angularAcceptance = []
for cos_bin in cos_bins:
    if cos_bin <= 0.:
        angularAcceptance.append(0.)
    else:
        angularAcceptance.append(cos_bin * winstonRatioFunc(-cos_bin))
angularAcceptance = clsim.I3CLSimFunctionFromTable(cos_bins[0], cos_bins[1]-cos_bins[0], angularAcceptance)




# extra DOM simulation for MICA DOMs (these were ignored in the tray segment)
tray.AddModule("I3PhotonToMCHitConverterForMDOMs", "mDOM_PMTs",
               RandomService = randomService,
               InputPhotonSeriesMapName = "PhotonsOnDOMs",
               OutputMCHitSeriesMapName = "MCHitSeriesMap_MICA_mDOM", 
               MCTreeName = "I3MCTree",
               IgnoreSubdetectors = ["AMANDA", "IceCube", "DeepCore", "IceTop", "PINGU"],
               PMTWavelengthAcceptance = wavelengthAcceptance,
               PMTAngularAcceptance = angularAcceptance,
               GlassAbsorptionLength = glassAbsLen,
               GlassThickness = glassThickness,
               GelAbsorptionLength = gelAbsLen
               )


# domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = 0.16510*I3Units.m)
# domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=True)
# pmtPhotonSimulator = clsim.I3CLSimPMTPhotonSimulatorIceCube(jitter=2.*I3Units.ns,
#                                                             pre_pulse_probability=0.,
#                                                             late_pulse_probability=0.,
#                                                             after_pulse_probability=0.)
# tray.AddModule("I3PhotonToMCHitConverter", "make_hits_IceCubDOM",
#                RandomService = randomService,
#                MCTreeName = "I3MCTree",
#                InputPhotonSeriesMapName = "PhotonsOnDOMs",
#                OutputMCHitSeriesMapName = "MCHitSeriesMap_MICA_stdIceCubeDOM",
#                #MCTreeName = "I3MCTree",
#                DOMRadiusWithoutOversize=0.16510*I3Units.m,
#                DOMOversizeFactor = 1.,
#                WavelengthAcceptance = domAcceptance,
#                AngularAcceptance = domAngularSensitivity,
#                PMTPhotonSimulator = pmtPhotonSimulator, # simulate jitter, no after-pulses and late-pulses
#                IgnoreDOMsWithoutDetectorStatusEntry = False,
#                DefaultRelativeDOMEfficiency=1.,
#                ReplaceRelativeDOMEfficiencyWithDefault=True)

# write the results to disk
tray.AddModule("I3Writer", "writer",
    filename = options.OUTFILE)



tray.Execute(options.NUMEVENTS+3)

