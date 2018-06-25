#!/usr/bin/env python
"""
A simple module to apply IC86 calibration and feature
extraction. This will also retrieve a bad DOM list from
the database and clean their pulses if there should be any.
"""

from __future__ import print_function
from optparse import OptionParser
import os
import string

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default=None,
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default="test_flashes_clsim_detsim.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("--simdata", action="store_true", default=False,
                  dest="SIMDATA", help="Use special wavedeform mode for DOMsimulator")
parser.add_option("--no-dynamic-bad-doms", action="store_false", default=True,
                  dest="DYNAMICBADDOMS", help="Do not retrieve a dynamic bad dom list from the database")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

########################
if options.INFILE:
        filename = options.INFILE
        if os.access(filename,os.R_OK) == False:
                raise RuntimeError("cannot find input file!")
        infile = filename
        print('using input file %s' % infile)
else:
        print("No input file!")
        parser.print_help()
        exit(-1)

infileRoot, infileExt = os.path.splitext(infile)
if infileExt == ".gz":
    infileRoot2, infileExt2 = os.path.splitext(infileRoot)
    if infileExt2 == ".i3":
        infileRoot=infileRoot2
        infileExt = ".i3.gz"

if infileExt != ".i3" and infileExt != ".i3.gz":
        raise Exception("you have to specify either a .i3 or an .i3.gz file!")

########################
outdir=""
outfile=None
if options.OUTFILE:
        outfile = options.OUTFILE
        # did the user specify a directory? then use that and auto-generate
        if os.path.isdir(outfile):
            outdir = outfile
            outfile = None
        else:
            outdir, outfile = os.path.split(outfile)

# add a trailing slash to the output directory name if not already there
if outdir and outdir!="":
    if outdir[-1] != "/":
        outdir += "/"

if not outfile:
        # automatically generate the output filename
        infileRootDir, infileRootFile = os.path.split(infileRoot)
        outfile = infileRootFile + "_calib"
        outfile = outfile + infileExt
print("output dir is %s" % outdir)
print("output file is %s" % outdir + outfile)

########################


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services

from icecube.sim_services import bad_dom_list_static
from icecube import WaveCalibrator
#from icecube import DomTools
icetray.load("DomTools", False)
from icecube import wavedeform


tray = I3Tray()

tray.AddModule("I3Reader","reader",
               Filename=infile)

if options.DYNAMICBADDOMS:
    dbserver="dbs2.icecube.wisc.edu"
    username="www"

    from icecube.sim_services.bad_dom_list_static import IC86_static_bad_dom_list, IC86_static_bad_dom_list_HLC
    xmlfile = os.path.expandvars('$I3_BUILD') + '/BadDomList/resources/scripts/QueryConfiguration.xml'
    staticBadDomListOnly = False

    tray.AddService( 'I3BadDomListFactory', 'BadDomListServiceHLC',
                     ServiceName = 'BadDomListServiceHLC',
                     Hostname = dbserver,
                     Timeout = 180,
                     QFileName = xmlfile
                     )
    tray.AddService( 'I3BadDomListFactory', 'BadDomListServiceSLC',
                     ServiceName = 'BadDomListServiceSLC',
                     Hostname = dbserver,
                     Timeout = 180,
                     QFileName = xmlfile
                     )

    badOms = IC86_static_bad_dom_list_HLC()
    tray.AddModule( 'I3BadDomListModule', 'BadDomsHLC',
                    BadDomListServiceName = 'BadDomListServiceHLC',
                    CleanedKeys = badOms,
                    BadDomsListVector = 'BadDomsList',
                    DbOverridesStatic = False,
                    StaticOverridesDb = staticBadDomListOnly
                    )

    badOmsSLC =IC86_static_bad_dom_list()
    tray.AddModule( 'I3BadDomListModule', 'BadDomsSLC',
                    BadDomListServiceName = 'BadDomListServiceSLC',
                    CleanedKeys = badOmsSLC,
                    BadDomsListVector = 'BadDomsListSLC',
                    DbOverridesStatic = False,
                    StaticOverridesDb = staticBadDomListOnly
                    )

if options.SIMDATA:
    bad_doms = bad_dom_list_static.IC86_static_bad_dom_list()
else:
    bad_doms = [OMKey(38,59), # Blackberry
                OMKey(68,42)  # Krabba
               ]

# online cleaning
tray.AddModule("I3DOMLaunchCleaning", 'baddomclean',
           CleanedKeys=bad_doms,
           InIceOutput='CleanInIceRawData',
           IceTopOutput='CleanIceTopRawData'
           )

# offline cleaning
tray.AddModule( 'I3DOMLaunchCleaning', 'OfflineLaunchCleaning',
    InIceInput = 'CleanInIceRawData',
    IceTopInput = 'CleanIceTopRawData',
    InIceOutput = 'OfflineCleanInIceRawData',
    IceTopOutput = 'OfflineCleanIceTopRawData',
    FirstLaunchCleaning = False,
    CleanedKeysList = 'BadDomsList'
    )


if not options.SIMDATA:
    print('Wavecalibrator for data selected')
    tray.AddModule('I3WaveCalibrator', 'wavecal',
               Launches='OfflineCleanInIceRawData',
               )
    tray.AddModule('I3Wavedeform', 'wavedeform',
           Output='UncleanedInIcePulses')
else:
    print('Wavecalibrator for DOMsimulator selected')
    tray.AddSegment(WaveCalibrator.DOMSimulatorCalibrator, 'wavecal',
               Launches='OfflineCleanInIceRawData',
               )
    tray.AddModule('I3Wavedeform', 'wavedeform',
           Output='UncleanedInIcePulses', UseDOMsimulatorTemplates=True)

# cleanup
tray.AddModule('Delete', 'finalcleanup',
    Keys = ["CalibratedWaveforms"])

# now we generate P-frames
tray.AddModule("I3NullSplitter", "nullSplit")

tray.AddModule("I3LCPulseCleaning", "I3LCPulseCleaning",
    Input='UncleanedInIcePulses',
    OutputHLC='UncleanedInIcePulsesHLC',
    OutputSLC='UncleanedInIcePulsesSLC')


tray.AddModule("I3Writer","writer",
    Filename = outdir+outfile)



tray.Execute()

