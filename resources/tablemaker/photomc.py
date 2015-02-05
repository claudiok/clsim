
from optparse import OptionParser
from icecube.icetray import I3Units
from os import path, unlink

usage = "usage: %prog [options] outputfile"
parser = OptionParser(usage)

parser.add_option("--seed", dest="seed", type="int", default=None,
    help="Seed for random number generators; harvested from /dev/random if unspecified.")
parser.add_option("--nevents", dest="nevents", type="int", default=100,
    help="Number of light sources to inject [%default]")
parser.add_option("--z", dest="z", type="float", default=0.,
    help="IceCube z-coordinate of light source, in meters [%default]")
parser.add_option("--zenith", dest="zenith", type="float", default=0.,
    help="Zenith angle of source, in IceCube convention and degrees [%default]")
parser.add_option("--energy", dest="energy", type="float", default=1,
    help="Energy of light source, in GeV [%default]")
parser.add_option("--light-source", choices=('cascade', 'flasher', 'infinite-muon'), default='cascade',
    help="Type of light source. If 'infinite-muon', Z will be ignored, and tracks sampled over all depths. [%default]")
parser.add_option("--step", dest="steplength", type="float", default=1,
    help="Sampling step length in meters [%default]")
parser.add_option("--overwrite", dest="overwrite", action="store_true", default=False,
    help="Overwrite output file if it already exists")
    
opts, args = parser.parse_args()

if len(args) != 1:
	parser.error("You must specify an output file!")
outfile = args[0]
if path.exists(outfile):
	if opts.overwrite:
		unlink(outfile)
	else:
		parser.error("Output file exists! Pass --overwrite to overwrite it.")

from icecube import icetray
from icecube.clsim.tablemaker.tabulator import CombinedPhotonGenerator, generate_seed

outfile = args[0]
if opts.seed is None:
	opts.seed = generate_seed()
opts.zenith *= I3Units.degree

from I3Tray import I3Tray

tray = I3Tray()

icetray.logging.set_level_for_unit('I3CLSimStepToTableConverter', 'TRACE')
icetray.logging.set_level_for_unit('I3CLSimTabulatorModule', 'DEBUG')
icetray.logging.set_level_for_unit('I3CLSimLightSourceToStepConverterGeant4', 'TRACE')
icetray.logging.set_level_for_unit('I3CLSimLightSourceToStepConverterFlasher', 'TRACE')

tray.AddSegment(CombinedPhotonGenerator, 'generator', Seed=opts.seed, PhotonSource=opts.light_source,
    Zenith=opts.zenith, ZCoordinate=opts.z, Energy=opts.energy, NEvents=opts.nevents, Filename=outfile)
    
tray.AddModule('TrashCan', 'MemoryHole')
tray.Execute()
tray.Finish()
