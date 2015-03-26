
all: timing_distributions_unweighted_mDOMs.pdf

# icemodels = lea lea_notilt lea_noanisotropy lea_notilt_noanisotropy
# ppc := $(patsubst %,test_events_ppc_%.i3,$(icemodels))
# clsim := $(patsubst %,test_events_clsim_%.i3,$(icemodels))

i3 := test_events_mcpes_unweighted.i3 test_events_mcpes_oversize_5.i3
hdf5 := $(addsuffix .hdf5,$(i3))

testingGCD.i3:
	./generateTestingGeometry.py --xpos=257 --ypos=212 --zpos=-399

test_events.i3: testingGCD.i3
	./generateTestEvents.py --energy 1e3 --numevents=10 --xpos=257 --ypos=212 --zpos=-399

test_events_mcpes_unweighted.i3: test_events.i3
	./makeHits.py --unweighted --oversize 1. test_events.i3 test_events_mcpes_unweighted.i3

test_events_mcpes_oversize_2.i3: test_events.i3
	./makeHits.py --oversize 2. test_events.i3 test_events_mcpes_oversize_2.i3

test_events_mcpes_oversize_5.i3: test_events.i3
	./makeHits.py --oversize 5. test_events.i3 test_events_mcpes_oversize_5.i3


# $(clsim): test_events.i3
# 	$(eval icemodel = $(basename $(subst test_events_clsim_,,$@)))
# 	./applyCLSim.py -p10 -i test_events.i3 --icemodel=test_ice_models/${icemodel}/ -o $@
#
# $(ppc): test_events.i3
# 	$(eval icemodel = $(basename $(subst test_events_ppc_,,$@)))
# 	./applyPPC.py -p10 -i test_events.i3 --icemodel=test_ice_models/${icemodel}/ -o $@

%.hdf5: %
	./extractData.py $<

timing_distributions_unweighted_mDOMs.pdf: $(hdf5)
	./generatePlots.py


