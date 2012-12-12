#!/bin/bash
set -e

echo "************** generating GCD.."
./generateTestingGeometry.py --xpos=257 --ypos=212 --zpos=-399

echo "************** generating events.."
./generateTestEvents.py --numevents=1000 --xpos=257 --ypos=212 --zpos=-399

echo "************** running CLSim.."
time ./applyCLSim.py -p10 -i test_events.i3 --icemodel=test_ice_models/lea/                     -o test_events_clsim_lea.i3
time ./applyCLSim.py -p10 -i test_events.i3 --icemodel=test_ice_models/lea_noanisotropy/        -o test_events_clsim_lea_noanisotropy.i3
time ./applyCLSim.py -p10 -i test_events.i3 --icemodel=test_ice_models/lea_notilt/              -o test_events_clsim_lea_notilt.i3
time ./applyCLSim.py -p10 -i test_events.i3 --icemodel=test_ice_models/lea_notilt_noanisotropy/ -o test_events_clsim_lea_notilt_noanisotropy.i3

echo "************** running PPC.."
time ./applyPPC.py -i test_events.i3 --icemodel=test_ice_models/lea/                     -o test_events_ppc_lea.i3
time ./applyPPC.py -i test_events.i3 --icemodel=test_ice_models/lea_noanisotropy/        -o test_events_ppc_lea_noanisotropy.i3
time ./applyPPC.py -i test_events.i3 --icemodel=test_ice_models/lea_notilt/              -o test_events_ppc_lea_notilt.i3
time ./applyPPC.py -i test_events.i3 --icemodel=test_ice_models/lea_notilt_noanisotropy/ -o test_events_ppc_lea_notilt_noanisotropy.i3

echo "************** rextracting data.."
./extractData.py test_events_clsim_lea.i3
./extractData.py test_events_clsim_lea_noanisotropy.i3
./extractData.py test_events_clsim_lea_notilt.i3
./extractData.py test_events_clsim_lea_notilt_noanisotropy.i3

./extractData.py test_events_ppc_lea.i3
./extractData.py test_events_ppc_lea_noanisotropy.i3
./extractData.py test_events_ppc_lea_notilt.i3
./extractData.py test_events_ppc_lea_notilt_noanisotropy.i3

echo "************** plot data.."

./generatePlots.py

echo "************** done!"


