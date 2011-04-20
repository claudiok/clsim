#!/bin/bash
set -e

echo "************** generating muons.."
time $I3_SRC/clsim/resources/scripts/generateTestMuons.py -o test_muons.i3 --apply-mmc --gcd=$I3_SRC/clsim/resources/GeoCalibDetectorStatus_IC86.55040_official.i3.gz -n1000

echo "************** running PPC.."
time $I3_SRC/ppc/resources/applyPPC.py -i test_muons.i3 -o test_muons_ppc.i3

echo "************** running CLSim.."
time $I3_SRC/clsim/resources/scripts/applyCLSim.py -i test_muons_ppc.i3 -o test_muons_ppc_clsim.i3 -m5 -p100

echo "************** booking.."
$I3_SRC/clsim/resources/scripts/bookThings.py test_muons_ppc_clsim.i3

echo "************** plotting.."
$I3_SRC/clsim/resources/scripts/plot_dom_occupancy.py test_muons_ppc_clsim.i3.hdf5

echo "************** done!"


