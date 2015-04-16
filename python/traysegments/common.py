from __future__ import print_function

"""
Convenience functions for configuring CLSim components.
"""

def configureOpenCLDevices(UseGPUs=True, UseCPUs=False, OverrideApproximateNumberOfWorkItems=None, DoNotParallelize=True, UseOnlyDeviceNumber=None):
    # get OpenCL devices
    from icecube.clsim import I3CLSimOpenCLDevice
    from icecube.icetray import logging
    import string
    
    openCLDevicesRaw = [device for device in I3CLSimOpenCLDevice.GetAllDevices() if (device.gpu and UseGPUs) or (device.cpu and UseCPUs)]
    openCLDevices = []
    
    # (auto-)configure OpenCL devices
    for i, device in enumerate(openCLDevicesRaw):
        if UseOnlyDeviceNumber is not None and i != UseOnlyDeviceNumber:
            # skip all devices except for the selected one (if there is a selection)
            continue

        if string.count(device.device, 'Tesla') > 0 or string.count(device.device, 'GTX') > 0:
            # assume these are "fast", all others are "slow"
            device.useNativeMath=True
            if string.count(device.device, 'Tesla') > 0 or string.count(device.device, '580') > 0 or string.count(device.device, '680') > 0 or string.count(device.device, '980') > 0:
                # these cards should have enough ram to support this
                device.approximateNumberOfWorkItems=1024000
            else:
                device.approximateNumberOfWorkItems=102400
        elif string.count(device.device, 'Tahiti') > 0:
            device.useNativeMath=True
            device.approximateNumberOfWorkItems=102400*2
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
                    logging.log_warn("failed to split CPU device into individual cores %s %s [using full device with minimal number of work-items to (hopefully) disable parallelization]" % (device.platform, device.device), unit="clsim")
                    device.approximateNumberOfWorkItems=1
                    openCLDevices.append(device)
            except:
                logging.log_error("failed to split CPU device into individual cores %s %s [using full device with minimal number of work-items to (hopefully) disable parallelization]" % (device.platform, device.device), unit="clsim")
                device.approximateNumberOfWorkItems=1
                openCLDevices.append(device)
            
        else:
            openCLDevices.append(device)
    
    return openCLDevices

def parseIceModel(IceModelLocation, disableTilt=False):
    from os.path import exists, isdir, isfile, expandvars
    from icecube.clsim.MakeIceCubeMediumProperties import MakeIceCubeMediumProperties
    from icecube.clsim.MakeAntaresMediumProperties import MakeAntaresMediumProperties
    from icecube.clsim.MakeIceCubeMediumPropertiesPhotonics import MakeIceCubeMediumPropertiesPhotonics
    
    if IceModelLocation=="ANTARES":
        return MakeAntaresMediumProperties()
    
    if not exists(IceModelLocation):
        raise RuntimeError("The specified ice model path \"%s\" does not exist" % IceModelLocation)
    
    if isdir(IceModelLocation):
        # it's a PPC ice description directory
        mediumProperties = MakeIceCubeMediumProperties(iceDataDirectory=IceModelLocation, useTiltIfAvailable=not disableTilt)
    elif isfile(IceModelLocation):
        # it's a photonics ice description file
        mediumProperties = MakeIceCubeMediumPropertiesPhotonics(tableFile=IceModelLocation)
    else:
        raise RuntimeError("The specified ice model path \"%s\" is neither a directory nor a file." % IceModelLocation)
    
    return mediumProperties
