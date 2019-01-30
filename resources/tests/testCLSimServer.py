#!/usr/bin/env python

from icecube import clsim, icetray, dataclasses
import time
from multiprocessing import Process
from numpy.random import uniform
from numpy import testing

icetray.logging.set_level('INFO')

def dummy_photon(step):
    photon = clsim.I3CLSimPhoton()
    for attr in 'x', 'y', 'z', 'theta', 'phi', 'time', 'weight':
        setattr(photon, attr, getattr(step, attr))
    photon.numScatters = 3
    photon.omID = 52
    photon.stringID = 23
    return photon

def dummy_photon_history(photon):
    history = clsim.I3CLSimPhotonHistory()
    for i in range(photon.numScatters):
        history.append(dataclasses.I3Position(i,i+0.5,i+3.14), i)
    return history

class DummyConverter(clsim.I3CLSimStepToPhotonConverter):
    def __init__(self):
        super(DummyConverter, self).__init__()
        self.input_queue = list()
    def IsInitialized(self):
        return True
    def GetWorkgroupSize(self):
        return 1
    def GetMaxNumWorkitems(self):
        return 64
    def EnqueueSteps(self, steps, id):
        self.input_queue.append((steps, id))
    def GetConversionResult(self):
        steps, id = self.input_queue.pop(0)
        photons = clsim.I3CLSimPhotonSeries(map(dummy_photon, steps))
        history = clsim.I3CLSimPhotonHistorySeries(map(dummy_photon_history, photons))
        return clsim.I3CLSimStepToPhotonConverter.ConversionResult_t(id, photons, history)

def test_client(client, num_steps):
    
    input_steps = []
    for i in range(num_steps):
        steps = clsim.I3CLSimStepSeries()
        for _ in range(int(uniform(10))+1):
            steps.append(clsim.I3CLSimStep())
            steps[-1].pos = dataclasses.I3Position(*uniform(size=3))
            steps[-1].dir = dataclasses.I3Direction(*uniform(size=3))
            steps[-1].time = uniform()
            steps[-1].weight = uniform()
        input_steps.append(steps)
        client.EnqueueSteps(steps, i)
    for i in range(num_steps):
        result = client.GetConversionResult()
        input = input_steps[result.identifier]
        assert len(result.photons) == len(input)
        assert len(result.photonHistories) == len(input)
        for step, photon, history in zip(input, result.photons, result.photonHistories):
        
            testing.assert_equal( photon.numScatters, 3 )
            testing.assert_equal( photon.omID, 52 )
            testing.assert_equal( photon.stringID, 23 )
            for attr in 'x', 'y', 'z', 'theta', 'phi', 'time', 'weight':
                testing.assert_equal( getattr(step, attr), getattr(photon, attr), err_msg='{} not equal'.format(attr))
            
            dummy_history = dummy_photon_history(photon)
            testing.assert_equal( len(dummy_history), len(history) )
            for pos, expected_pos in zip(history, dummy_history):
                testing.assert_equal( pos, expected_pos )

# First, ensure that the test passes when the converter is called directly
test_client(DummyConverter(), 10)

# Now, call through the server in a separate process
converters = clsim.I3CLSimStepToPhotonConverterSeries([DummyConverter()])
address = 'ipc:///tmp/clsim-server.ipc'
server = clsim.I3CLSimServer(address,converters)

def fire_a_few():
    client = clsim.I3CLSimClient(address)    
    testing.assert_equal( client.workgroupSize, 1 )
    testing.assert_equal( client.maxNumWorkitems, 64 )
    test_client(client, 10)

procs = [Process(target=fire_a_few) for i in range(3)]
for p in procs:
    p.start()
for p in procs:
    p.join()
    assert p.exitcode == 0

print("going to destroy")
del server
print("destroyed")


