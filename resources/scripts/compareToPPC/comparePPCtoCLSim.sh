#!/bin/bash
set -e

echo "************** generating muons.."
time $I3_BUILD/clsim/resources/scripts/compareToPPC/generateTestMuons.py -o test_muons.i3 --gcd=$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected.i3.gz -n1000

echo "************** running PPC.."
time $I3_BUILD/clsim/resources/scripts/compareToPPC/applyPPC.py -i test_muons.i3 -o test_muons_ppc.i3

echo "************** running CLSim.."
time $I3_BUILD/clsim/resources/scripts/compareToPPC/applyCLSim.py -i test_muons_ppc.i3 -o test_muons_ppc_clsim.i3 -p100

echo "************** booking.."
$I3_BUILD/clsim/resources/scripts/compareToPPC/bookThings.py test_muons_ppc_clsim.i3

echo "************** plotting.."
$I3_BUILD/clsim/resources/scripts/compareToPPC/plot_dom_occupancy.py test_muons_ppc_clsim.i3.hdf5

echo "************** done!"


