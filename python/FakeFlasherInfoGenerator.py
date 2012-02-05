from .. import icetray, dataclasses

from I3Tray import I3Units

class FakeFlasherInfoGenerator(icetray.I3ConditionalModule):
    """
    Generate fake I3FlasherInfoVect frame objects.
    """
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("FlasherInfoVectName",
                          "Name of the I3FlasherInfoVect object to write",
                          "I3FlasherInfo")
        self.AddParameter("FlashingDOM",
                          "OMKey of the flashing DOM",
                          icetray.OMKey(0,0))
        self.AddParameter("FlasherTime",
                          "Time of flash, relative to event start time",
                          0.*I3Units.ns)
        self.AddParameter("FlasherMask",
                          "bitmask (12bits) describing the active flashers (index 0-11)",
                          4032) # Only the 6 horizontal LEDs (111111000000)
        self.AddParameter("FlasherBrightness",
                          "flasher brightness setting (0-127)",
                          127)
        self.AddParameter("FlasherWidth",
                          "flasher width setting (0-127)",
                          127)
                          
        self.AddOutBox("OutBox")
        
    def Configure(self):
        self.flasherInfoVectName = self.GetParameter("FlasherInfoVectName")
        self.flashingDOM = self.GetParameter("FlashingDOM")
        self.flasherTime = self.GetParameter("FlasherTime")
        self.flasherMask = self.GetParameter("FlasherMask")
        self.flasherBrightness = self.GetParameter("FlasherBrightness")
        self.flasherWidth = self.GetParameter("FlasherWidth")
    
    def DAQ(self, frame):
        # generate the output vector
        outputVect = dataclasses.I3FlasherInfoVect()
        
        # fill an element
        flasherInfo = dataclasses.I3FlasherInfo()
        
        # compatible to "photoflash"
        flasherInfo.flashing_om = self.flashingDOM
        flasherInfo.led_brightness = self.flasherBrightness
        flasherInfo.mask = self.flasherMask
        flasherInfo.width = self.flasherWidth
        flasherInfo.flash_time = self.flasherTime
        
        # Flasher info default (dummy) values
        flasherInfo.rate = 0
        flasherInfo.atwd_bin_size = 0
        
        # add the element
        outputVect.append(flasherInfo)

        # add to frame and push
        frame[self.flasherInfoVectName] = outputVect
        self.PushFrame(frame)
