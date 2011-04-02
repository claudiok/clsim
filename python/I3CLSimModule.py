from icecube import icetray, dataclasses

from icecube.clsim import I3CLSimStepToPhotonConverterOpenCL
from icecube.clsim import I3CLSimParticleToStepConverterGeant4
from icecube.clsim import I3CLSimParticleToStepConverterGeant4MP
from icecube.clsim import I3CLSimSimpleGeometryFromI3Geometry
from icecube.clsim import I3CLSimPhotonSeries, I3CLSimStepSeries
from icecube.clsim import I3CLSimPhotonSeriesToPhotonSeriesMapConverter

import string

from I3Tray import I3Units

def initializeOpenCL(rng, geometry, medium):
    conv = I3CLSimStepToPhotonConverterOpenCL(RandomService=rng,
                                                    UseNativeMath=True)
    
    print "available OpenCL devices:"
    deviceList = conv.GetDeviceList()
    for device in deviceList:
        print "platform: \"%s\", device: \"%s\"" % (device[0], device[1])
        
    print ""
    
    # do a semi-smart device selection
    deviceToUse=None
    
    # look for a "GeForce" device first
    geForceDevices = [device for device in deviceList if string.lower(device[1]).find("geforce")!=-1]
    if len(geForceDevices)>0:
        if len(geForceDevices)>1:
            deviceToUse=geForceDevices[0]
            print "You seem to have more than one GeForce GPU. Using the first one. (\"%s\")" % (deviceToUse[1])
        else:
            deviceToUse=geForceDevices[0]
            print "You seem to have a GeForce GPU. (\"%s\")" % (deviceToUse[1])
    else:
        print "No GeForce device found. Just using the first available one."
        deviceToUse=deviceList[0]
    
    print " -> using OpenCL device \"%s\" on platform \"%s\"." % (deviceToUse[1], deviceToUse[0])
    
    conv.SetDeviceName(deviceToUse[0], deviceToUse[1])

    conv.SetMediumProperties(medium)
    conv.SetGeometry(geometry)

    conv.Compile()
    #print conv.GetFullSource()

    conv.workgroupSize = conv.maxWorkgroupSize
    
    # use approximately 512000 work items, convert to a multiple of the workgroup size
    conv.maxNumWorkitems = (512000/conv.workgroupSize)*conv.workgroupSize

    print "maximum workgroup size is", conv.maxWorkgroupSize
    print "configured workgroup size is", conv.workgroupSize
    print "maximum number of work items is", conv.maxNumWorkitems

    conv.Initialize()

    return conv
    
def initializeGeant4(rng, medium, openCLconverter, multiprocessor=False):
    if not multiprocessor:
        conv = I3CLSimParticleToStepConverterGeant4(rng.Integer(900000000))
    else:
        randomSeeds=[]
        # we will NOT have more than 100 cores.. ;-)
        for i in range(100): randomSeeds.append(rng.Integer(900000000))
        conv = I3CLSimParticleToStepConverterGeant4MP(randomSeeds=randomSeeds, numInstances=4)

    conv.SetMediumProperties(medium)
    conv.SetMaxBunchSize(openCLconverter.maxNumWorkitems)
    conv.SetBunchSizeGranularity(openCLconverter.workgroupSize)
    
    #conv.SetElectronPositronMinEnergyForSecondary(1.*I3Units.GeV)
    #conv.SetElectronPositronMaxEnergyForSecondary(10.*I3Units.TeV)
    
    conv.Initialize()
    
    return conv



class I3CLSimModule(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('RandomService',
                          'A random number generating service (derived from I3RandomService).',
                          None)
        self.AddParameter('MediumProperties',
                          'An instance of I3CLSimMediumProperties describing the ice/water proerties.',
                          None)
        self.AddParameter('MaxNumParallelEvents',
                          'Maximum number of events that will be processed by the GPU in parallel.',
                          1000)
        self.AddParameter('MCTreeName',
                          'Name of the I3MCTree frame object. All particles except neutrinos will be read from this tree.',
                          "I3MCTree")
        self.AddParameter('PhotonSeriesMapName',
                          'Name of the I3CLSimPhotonSeriesMap frame object that will be written to the frame.',
                          "I3CLSimPhotonSeriesMap")
        self.AddParameter('IgnoreMuons',
                          'If set to True, muons will not be propagated.',
                          False)
        
        self.neutrinoTypes = [dataclasses.I3Particle.NuE,
                              dataclasses.I3Particle.NuEBar,
                              dataclasses.I3Particle.NuMu,
                              dataclasses.I3Particle.NuMuBar,
                              dataclasses.I3Particle.NuTau,
                              dataclasses.I3Particle.NuTauBar]

        self.muonTypes = [dataclasses.I3Particle.MuMinus,
                          dataclasses.I3Particle.MuPlus]

        self.frameCache = dict()
    
    def Configure(self):
        print "Entering I3CLSimModule::Configure()"
        
        self.randomService = self.GetParameter('RandomService')
        self.mediumProperties = self.GetParameter('MediumProperties')
        self.maxNumParallelEvents = self.GetParameter('MaxNumParallelEvents')
        self.MCTreeName = self.GetParameter('MCTreeName')
        self.photonSeriesMapName = self.GetParameter('PhotonSeriesMapName')
        self.ignoreMuons = self.GetParameter('IgnoreMuons')

        if self.randomService is None: raise RuntimeError("You have to specify the \"RandomService\" parameter!")
        if self.mediumProperties is None: raise RuntimeError("You have to specify the \"MediumProperties\" parameter!")
        if self.maxNumParallelEvents <= 0: raise RuntimeError("Values <= 0 are invalid for the \"MaxNumParallelEvents\" parameter.!")
        
        self.currentEventID = 0
        self.geometryIsConfigured = False
        self.totalSimulatedEnergyForFlush = 0.
        self.totalNumParticlesForFlush = 0
    
        self.photonSeriesMapConverter = I3CLSimPhotonSeriesToPhotonSeriesMapConverter()
    
    def Geometry(self, frame):
        print "Entering I3CLSimModule::Geometry()"
        
        if self.geometryIsConfigured:
            raise RuntimeError("This module currently supports only a single geometry per input file.")
        
        print "Retrieving geometry.."
        self.I3Geometry = frame["I3Geometry"]
                
        print "Converting geometry.."
        self.geometry = I3CLSimSimpleGeometryFromI3Geometry(43.18*I3Units.cm/2., self.I3Geometry)
                
        print "Initializing CLSim.."
        # initialize OpenCL
        self.openCLStepsToPhotonsConverter = initializeOpenCL(self.randomService,
                                                              self.geometry,
                                                              self.mediumProperties)

        print "Initializing Geant4.."
        # initialize Geant4 (will set bunch sizes according to the OpenCL settings)
        self.geant4ParticleToStepsConverter = initializeGeant4(self.randomService,
                                                               self.mediumProperties,
                                                               self.openCLStepsToPhotonsConverter,
                                                               multiprocessor=False) # the multiprocessor version is not yet safe to use

        print "Initializing photon map converter.."
        self.photonSeriesMapConverter.SetGeometry(self.geometry)
                
        print "Initialization complete."
        
        self.geometryIsConfigured=True
        self.PushFrame(frame)
    
    def _FlushFrameCache(self):
        print "Flushing frame cache.."

        self.geant4ParticleToStepsConverter.EnqueueBarrier()
        
        numBunchesSentToOpenCL = 0
        #while self.geant4ParticleToStepsConverter.BarrierActive() or self.geant4ParticleToStepsConverter.MoreStepsAvailable():
        while True:
            if not self.geant4ParticleToStepsConverter.MoreStepsAvailable():
                #print "No steps are available right now"
                if not self.geant4ParticleToStepsConverter.BarrierActive():
                    #print "No barrier is active, we are done"
                    break # nothing to read and no barrier. We are done
            
                # no steps available, but barrier is still enqueued. Geant4
                # must still be working. Try to get steps, but timeout after a while.
                #print "Waiting until steps become available or a timeout occurs"
                steps = self.geant4ParticleToStepsConverter.GetConversionResult(0.1*I3Units.s) # time in seconds
                if steps is None:
                    #print "Timeout while waiting for steps, trying again"
                    continue
                #print "Steps retrieved"
            else:
                # there are steps available, fetch them!
                #print "Steps available, retrieving them"
                steps = self.geant4ParticleToStepsConverter.GetConversionResult()
                #print "Steps retrieved"
            
            if isinstance(steps, tuple) and isinstance(steps[1], dataclasses.I3Particle):
                print "Got a secondary with E =", steps[1].energy, "(ignoring it: needs to be fixed!)"
                #print "Got a secondary particle from Geant4. This was not configured, something is wrong. ignoring."
                continue
            elif not isinstance(steps, I3CLSimStepSeries):
                print "Got something unknown from Geant4. ignoring."
                continue
            
            # has to be a bunch of steps at this point
            print "Got", len(steps), "steps from Geant4, sending them to OpenCL"
            
            self.openCLStepsToPhotonsConverter.EnqueueSteps(steps, numBunchesSentToOpenCL)
            numBunchesSentToOpenCL += 1

        print "Geant4 finished, retrieving results from GPU.."

        # join photon lists
        allPhotons = I3CLSimPhotonSeries()
        for i in range(numBunchesSentToOpenCL):
            #while openCLStepsToPhotonsConverter.MorePhotonsAvailable():
            res = self.openCLStepsToPhotonsConverter.GetConversionResult()
            photons = res[1]
            print "result #", res[0], ": num photons =", len(res[1])
            if (len(photons)>0):
                allPhotons.extend(photons) # append the photons to the list of all photons
            del photons
            del res

        print "results fetched from OpenCL."

        print "Got", len(allPhotons), "photons in total during flush."

        print "finished."
    
        self.photonSeriesMapConverter.SetPhotonSeries(allPhotons)

        for identifier, frame in sorted(self.frameCache.items()):
            print "pushing frame number", identifier, "..."
            
            photonSeriesMap = self.photonSeriesMapConverter.GetPhotonSeriesMapForIdentifier(identifier)
            
            frame[self.photonSeriesMapName] = photonSeriesMap
            
            self.PushFrame(frame)

        # TODO: add a sanity check: are there identifiers in the photon series that were not retrieved?
                
        self.frameCache = dict() # just empty the cache for now

    
    def Physics(self, frame):
        #print "Entering I3CLSimModule::Physics()"
        if not self.geometryIsConfigured:
            raise RuntimeError("Received Physics frame before Geometry frame")
        
        MCTree = frame[self.MCTreeName]
        
        for particle in MCTree:
            if particle.type in self.neutrinoTypes: continue
            if (self.ignoreMuons) and (particle.type in self.muonTypes): continue
            
            self.totalSimulatedEnergyForFlush += particle.energy
            self.totalNumParticlesForFlush += 1
            self.geant4ParticleToStepsConverter.EnqueueParticle(particle, self.currentEventID)

        self.frameCache[self.currentEventID] = frame
        #self.PushFrame(frame)
        
                
        if len(self.frameCache) >= self.maxNumParallelEvents:
            print "Flushing results for a total energy of", self.totalSimulatedEnergyForFlush/I3Units.GeV, "GeV for", self.totalNumParticlesForFlush, "particles"
            self.totalSimulatedEnergyForFlush=0.
            self.totalNumParticlesForFlush=0
            
            self._FlushFrameCache()
                
        # new event ID
        self.currentEventID += 1
        
    def Finish(self):
        print "Entering I3CLSimModule::Finish()"

        print "Flushing results for a total energy of", self.totalSimulatedEnergyForFlush/I3Units.GeV, "GeV for", self.totalNumParticlesForFlush, "particles"
        self.totalSimulatedEnergyForFlush=0.
        self.totalNumParticlesForFlush=0
        self._FlushFrameCache()

        print "Flushing I3Tray.."
        self.Flush()

