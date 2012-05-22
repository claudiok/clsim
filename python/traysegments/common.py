
"""
Convenience functions for configuring CLSim components.
"""

def configureOpenCLDevices(UseGPUs=True, UseCPUs=False, OverrideApproximateNumberOfWorkItems=None, DoNotParallelize=True):
    # get OpenCL devices
    from icecube.clsim import I3CLSimOpenCLDevice
    import string
    
    openCLDevicesRaw = [device for device in I3CLSimOpenCLDevice.GetAllDevices() if (device.gpu and UseGPUs) or (device.cpu and UseCPUs)]
    openCLDevices = []
    
    # (auto-)configure OpenCL devices
    for device in openCLDevicesRaw:
        if string.count(device.device, 'Tesla') > 0 or string.count(device.device, 'GTX') > 0:
            # assume these are "fast", all others are "slow"
            device.useNativeMath=True
            if string.count(device.device, 'Tesla') > 0 or string.count(device.device, '580') > 0 or string.count(device.device, '680') > 0:
                # these cards should have enough ram to support this
                device.approximateNumberOfWorkItems=1024000
            else:
                device.approximateNumberOfWorkItems=102400
        else:
            device.useNativeMath=False
            device.approximateNumberOfWorkItems=10240
            
        if OverrideApproximateNumberOfWorkItems is not None:
            device.approximateNumberOfWorkItems=OverrideApproximateNumberOfWorkItems

        if DoNotParallelize and device.cpu:
            # check if we can split this device into individual cores
            try:
                if device.platform == "Intel(R) OpenCL":
                    # device fission seems to cause serious segfaults in the Intel OpenCL driver.
                    # do not use it.
                    subDevices = []
                else:
                    subDevices = device.SplitDevice()
                
                if len(subDevices) > 0:
                    if OverrideApproximateNumberOfWorkItems is not None:
                        subDevices[0].approximateNumberOfWorkItems=OverrideApproximateNumberOfWorkItems
                    else:
                        subDevices[0].approximateNumberOfWorkItems=10240
                    openCLDevices.append(subDevices[0])
                else:
                    print "failed to split CPU device into individual cores", device.platform, device.device, "[using full device with minimal number of work-items to (hopefully) disable parallelization]"
                    device.approximateNumberOfWorkItems=1
                    openCLDevices.append(device)
            except:
                print "failed to split CPU device into individual cores", device.platform, device.device, "(exception) [using full device with minimal number of work-items to (hopefully) disable parallelization]"
                device.approximateNumberOfWorkItems=1
                openCLDevices.append(device)
            
        else:
            openCLDevices.append(device)
    
    return openCLDevices

def parseIceModel(IceModelLocation):
    from os.path import exists, isdir, expandvars
    from icecube.clsim.MakeIceCubeMediumProperties import MakeIceCubeMediumProperties
    from icecube.clsim.MakeIceCubeMediumPropertiesPhotonics import MakeIceCubeMediumPropertiesPhotonics
    
    if not exists(IceModelLocation):
        raise RuntimeError("The specified ice model path \"%s\" does not exist" % IceModelLocation)
    
    if isdir(IceModelLocation):
        # it's a PPC ice description directory
        mediumProperties = MakeIceCubeMediumProperties(iceDataDirectory=IceModelLocation)
    elif isfile(IceModelLocation):
        # it's a photonics ice description file
        mediumProperties = MakeIceCubeMediumPropertiesPhotonics(tableFile=IceModelLocation)
    else:
        raise RuntimeError("The specified ice model path \"%s\" is neither a directory nor a file." % IceModelLocation)
    
    return mediumProperties
