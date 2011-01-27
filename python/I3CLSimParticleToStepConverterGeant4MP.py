from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimParticleToStepConverter, I3CLSimParticleToStepConverterGeant4, I3CLSimMediumProperties
from icecube.clsim import I3CLSimStepSeries

### TODO: seed the rngs!

import multiprocessing as mp
from Queue import Empty

def checkGeant4Result(result, prefix=""):
    if (prefix is not None) and (prefix != ""):
        print prefix,
    if isinstance(result, tuple) and isinstance(result[1], dataclasses.I3Particle):
        if prefix is not None:
            print "Got a particle from Geant4 with id", result[0], "and energy", result[1].GetEnergy()/I3Units.GeV, "GeV"
        return True
    elif isinstance(result, I3CLSimStepSeries):
        if prefix is not None:
            print "Got", len(result), "steps from Geant4"
        return True
    else:
        if prefix is not None:
            print "Got something unknown from Geant4."
        return False


def TheGeant4Process(procNum, 
                     queueFromParent,
                     queueToParent,
                     maxQueueItems,
                     randomSeed,
                     physicsListName,
                     maxBetaChangePerStep,
                     maxNumPhotonsPerStep,
                     medium,
                     maxBunchSize,
                     bunchSizeGranularity,
                     minESecondary,
                     maxESecondary,
                     noisy):
    if noisy: print "worker #", procNum, "started."
    
    # set up Geant4
    conv = I3CLSimParticleToStepConverterGeant4(maxQueueItems)
    conv.SetMediumProperties(medium)
    conv.SetMaxBunchSize(maxBunchSize)
    conv.SetBunchSizeGranularity(bunchSizeGranularity)
    conv.electronPositronMinEnergyForSecondary = minESecondary
    conv.electronPositronMaxEnergyForSecondary = maxESecondary
    conv.Initialize()
    
    def checkWorkItem(workItem):
        if workItem is None: return True # None is ok, will be treated as barrier
        if not isinstance(workItem, tuple):
            #print "Got unknown work. Ignoring."
            return False
        if len(workItem) != 2:
            #print "Got a tuple of length", len(workItem), ", expected length 2. Ignoring."
            return False
        if not isinstance(workItem[0], dataclasses.I3Particle):
            #print "workItem[0] is not an instance of dataclasses.I3Particle. Ignoring."
            return False
        if not isinstance(workItem[1], int):
            #print "workItem[1] is not an int. Ignoring."
            return False
        return True
    
    if noisy:
        workerNamePrefixForLogging = "worker #%u" % procNum
    else:
        workerNamePrefixForLogging = None
    
    shouldFinish=False
    while True:
        workItem = queueFromParent.get()
        
        if isinstance(workItem, str) and workItem=="Terminate":
            if noisy: print "worker #", procNum, "received termination signal."
            break
        
        if not checkWorkItem(workItem):
            if noisy: print "worker #", procNum, "ignored invalid work item."
            continue
        
        if workItem is not None:
            if noisy: print "worker #", procNum, "about to add work item #", workItem[1]
            conv.EnqueueParticle(workItem[0], workItem[1])
            if noisy: print "worker #", procNum, "added work item #", workItem[1]
            
            # check if there is something on the output queue
            while conv.MoreStepsAvailable():
                if noisy: print "worker #", procNum, "seems to have something ready on the Geant4 result queue."
                s = conv.GetConversionResult()
                if not checkGeant4Result(s, workerNamePrefixForLogging):
                    print "worker #", procNum, "ignored invalid Geant4 result."
                    continue
                if noisy: print "worker #", procNum, "got result from Geant4, sending to parent."
                queueToParent.put(s)
                if noisy: print "worker #", procNum, "sent Geant4 result to parent."
                del s
            
            continue
        
        # if we got here, a barrier was received.
        # wait for all current Geant4 jobs to finish and add them
        # to the output queue.
        if noisy: print "worker #", procNum, "received barrier. waiting for Geant4 to finish processing.."
        conv.EnqueueBarrier()
        
        while conv.BarrierActive() or conv.MoreStepsAvailable():
            if noisy: print "worker #", procNum, "More steps available, receiving.."
            s = conv.GetConversionResult()
            if not checkGeant4Result(s, workerNamePrefixForLogging):
                print "worker #", procNum, "ignored invalid Geant4 result."
                continue
            if noisy: print "worker #", procNum, "got result from Geant4, sending to parent."
            
            secondItem=None
            if not conv.BarrierActive():
                if noisy: print "worker #", procNum, "finished processing barrier. adding notification along with the last result."
                secondItem=procNum
            
            queueToParent.put((s, secondItem))

            if noisy: print "worker #", procNum, "sent Geant4 result to parent."

            if conv.MoreStepsAvailable():
                continue

            # hack: apparently something in multiprocessing.Queue goes wrong when using two queues..
            # poking the _other_ queue seems to help the parent to receive the messages
            try:
                workItem = queueFromParent.get(True, 1.)
                if isinstance(workItem, str) and workItem=="Terminate":
                    if noisy: print "worker #", procNum, "received termination signal."
                    shouldFinish=True
                    break
                print "worker #", procNum, "received work from parent, but expected none (barrier is active). Ignoring."
            except Empty:
                pass
            
            
            
            #del s
        
        
        if shouldFinish: break
        if noisy: print "worker #", procNum, "finished processing barrier. re-starting receive loop."

    
    del conv
    
    if noisy: print "worker #", procNum, "finished."



class I3CLSimParticleToStepConverterGeant4MP(I3CLSimParticleToStepConverter):
    """
    A multiprocessing version of I3CLSimParticleToStepConverterGeant4.
    Uses forked processes, each one with their own I3CLSimParticleToStepConverterGeant4
    instance to parallelize Geant4 processing.
    """
    
    def _print_if_noisy(self, *args):
        if not self.noisy: return
        for arg in args:
            print arg,
        print ""
    
    def __init__(self, 
                 randomSeeds,
                 physicsListName="QGSP_BERT_EMV",
                 maxBetaChangePerStep=0.1,
                 maxNumPhotonsPerStep=200.,
                 maxQueueItems=20,
                 numInstances=mp.cpu_count(),
                 noisy=False):
        self.maxQueueItems = maxQueueItems
        
        self.randomSeeds=randomSeeds
        if isinstance(self.randomSeeds, int): self.randomSeeds=[self.randomSeeds]
        if len(self.randomSeeds) < numInstances:
            raise ValueError("You need to specify a random seed for each instance")
        
        self.physicsListName=physicsListName
        self.maxBetaChangePerStep=maxBetaChangePerStep
        self.maxNumPhotonsPerStep=maxNumPhotonsPerStep
        
        self.noisy=noisy
        
        self.queueFromProcesses = None
        self.numInstances = numInstances
        self.processes = []
        self.queuesToProcess = []
        self.currentProcess=0
        
        self.processHasBarrier = []
        
        self.initialized = False
        self.medium = None
        self.maxBunchSize = 51200
        self.bunchSizeGranularity = 512
        
        # set up a temporary instance of I3CLSimParticleToStepConverterGeant4 to get some default values
        tempInst = I3CLSimParticleToStepConverterGeant4(maxQueueItems)
        self.minESecondary = tempInst.electronPositronMinEnergyForSecondary
        self.maxESecondary = tempInst.electronPositronMaxEnergyForSecondary
        del tempInst
        
        self._print_if_noisy("construction complete.")
        
    def __del__(self):
        self._print_if_noisy("sending kill message to children..")
        for queue in self.queuesToProcess:
            queue.put("Terminate")
        self._print_if_noisy("waiting for children to terminate..")
        for i, process in enumerate(self.processes):
            if process.is_alive():
                process.join(5.) # 5 second timeout before we kill it anyway..
                if process.is_alive():
                    print " *** Process", i, "failed to terminate. Killing it. ***"
                    process.terminate()
        self._print_if_noisy("cleanup complete.")
        
    def SetElectronPositronMinEnergyForSecondary(self, val):
        if self.initialized: raise RuntimeError("Already initialized")
        self.minESecondary = val
        
    def SetElectronPositronMaxEnergyForSecondary(self, val):
        if self.initialized: raise RuntimeError("Already initialized")
        self.maxESecondary = val
    
    def GetElectronPositronMinEnergyForSecondary(self):
        if self.initialized: raise RuntimeError("Already initialized")
        return self.minESecondary

    def GetElectronPositronMaxEnergyForSecondary(self):
        if self.initialized: raise RuntimeError("Already initialized")
        return self.maxESecondary

    # make those properties, too
    electronPositronMinEnergyForSecondary = property(GetElectronPositronMinEnergyForSecondary, SetElectronPositronMinEnergyForSecondary)
    electronPositronMaxEnergyForSecondary = property(GetElectronPositronMaxEnergyForSecondary, SetElectronPositronMaxEnergyForSecondary)
    
    
    def SetBunchSizeGranularity(self, num):
        if self.initialized: raise RuntimeError("Already initialized")
        self.bunchSizeGranularity = num

    def SetMaxBunchSize(self, num):
        if self.initialized: raise RuntimeError("Already initialized")
        self.maxBunchSize = num

    def SetMediumProperties(self, medium):
        if self.initialized: raise RuntimeError("Already initialized")
        if not isinstance(medium, I3CLSimMediumProperties):
            raise ValueError("expected argument of type I3CLSimMediumProperties!")
        self.medium = medium
    
    def Initialize(self):
        if self.initialized: raise RuntimeError("Already initialized")
        if self.medium is None: raise RuntimeError("No medium configured. Use SetMediumProperties().")
        
        self._print_if_noisy("initializing..")
        self.queueFromProcesses = mp.Queue()
        
        # set up the processes
        self._print_if_noisy("setting up processes..")
        self.processes = []
        self.queuesToProcess = []
        self.processHasBarrier = []
        for i in range(self.numInstances):
            self.processHasBarrier.append(False)
            queueToChild = mp.Queue()
            self.queuesToProcess.append(queueToChild)
            self.processes.append( \
             mp.Process( \
              target=TheGeant4Process, \
              args=(i,\
                    queueToChild, \
                    self.queueFromProcesses, \
                    self.maxQueueItems, \
                    self.randomSeeds[i], \
                    self.physicsListName, \
                    self.maxBetaChangePerStep, \
                    self.maxNumPhotonsPerStep, \
                    self.medium, \
                    self.maxBunchSize, \
                    self.bunchSizeGranularity, \
                    self.minESecondary, \
                    self.maxESecondary, \
                    self.noisy, \
                   ) \
             ) \
            )
        
        self._print_if_noisy("starting processes..")
        for process in self.processes:
            process.daemon = True # terminate child if parent dies
            process.start()
        
        self._print_if_noisy("initialized!")
        self.initialized = True
        
    def IsInitialized(self):
        return self.initialized
        
    def EnqueueParticle(self, particle, identifier):
        if not self.initialized: raise RuntimeError("Not initialized!")
        if self.BarrierActive(): raise RuntimeError("Cannot enqeue particle while barrier is active")
        
        #self._print_if_noisy("about to send particle to process #", self.currentProcess)
        queue = self.queuesToProcess[self.currentProcess]
        
        queue.put((particle, identifier)) # put a tuple on the queue
        
        self.currentProcess += 1
        if self.currentProcess >= len(self.queuesToProcess):
            self.currentProcess = 0
        
    def EnqueueBarrier(self):
        if not self.initialized: raise RuntimeError("Not initialized!")
        if self.BarrierActive(): raise RuntimeError("Cannot enqeue barrier while another barrier is active")
        self._print_if_noisy("Enqueuing barrier..")
        for i in range(len(self.processHasBarrier)):
            self.processHasBarrier[i] = True
        for queue in self.queuesToProcess:
            queue.put(None)
        self._print_if_noisy("Barrier enqueued..")
        if not self.BarrierActive():
            raise RuntimeError("Internal error, barrier should be active")
        
    def BarrierActive(self):
        if not self.initialized: raise RuntimeError("Not initialized!")
        for val in self.processHasBarrier:
            if val: return True
        return False
        
    def MoreStepsAvailable(self):
        if not self.initialized: raise RuntimeError("Not initialized!")
        return not self.queueFromProcesses.empty()
        
    def GetConversionResult(self):
        if not self.initialized: raise RuntimeError("Not initialized!")
        
        self._print_if_noisy("Trying to get a conversion result..")
        
        if self.noisy:
            prefixForLogging = "parent"
        else:
            prefixForLogging = None
        
        while True:
            # get one item from the queue
            item = self.queueFromProcesses.get()
            #while True:
            #    try:
            #        item = self.queueFromProcesses.get(True, 1.)
            #        break
            #    except:
            #        print "queue size:", self.queueFromProcesses.qsize()
            #        pass
        
            self._print_if_noisy("Got something. queue size after receive:", self.queueFromProcesses.qsize())
        
            if not isinstance(item, tuple):
                print "Expected to receive a tuple from the child. Got something else. Ignoring."
                continue
            if len(item) != 2:
                print "Expected a tuple of length 2 from the child. Got length", len(item), ". Ignoring."
            if not item[1] is None:
                if not isinstance(item[1], int):
                    print "Expected the second tuple item to either be None or an integer. Got something else. Treating as None."
                else:
                    if item[1] >= len(self.processHasBarrier):
                        print "got barrier-finished notification from unknown process #", item[1], "Ignoring."
                        continue
                    if not self.processHasBarrier[item[1]]:
                        print "got barrier-finished notification process #", item[1], "which does not have a barrier enqueued. Ignoring."
                        continue
                    self._print_if_noisy("re-setting the barrier flag for process #", item[1])
                    self.processHasBarrier[item[1]] = False
            
            if not checkGeant4Result(item[0], prefixForLogging):
                print "Got something unknown from child. Ignoring."
                continue
            else:
                self._print_if_noisy("Got result from child. returning it now!")
                return item[0]



    
    
    