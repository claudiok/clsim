#
# Copyright (c) 2011, 2012
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 
# 
# $Id: FlasherInfoVectToFlasherPulseSeriesConverter.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file FlasherInfoVectToFlasherPulseSeriesConverter.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#


from .. import icetray, dataclasses
from . import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries

from I3Tray import I3Units
import math

class FlasherInfoVectToFlasherPulseSeriesConverter(icetray.I3ConditionalModule):
    """
    Read I3FlasherInfo objects from the frame, apply knowledge taken from
    various places (wiki, C.Wendt, ppc, photonics, ...) and create
    I3CLSimFLasherPulse objects desribing the light output of individual
    LEDs.
    """
    
    # Hard-coded list of color flashers (cDOMs) for IceCube 86.
    # All other DOMs are considered as "standard" DOMs.
    # (from http://wiki.icecube.wisc.edu/index.php/CDOM_Info)
    colorDOMs = set([icetray.OMKey(79,  1),  # UP9P3904  Frogs
                     icetray.OMKey(79,  8),  # TP9P3901  Goopy_Rain
                     icetray.OMKey(79, 13),  # UP9P3898  Torrential
                     icetray.OMKey(79, 22),  # TP9P3903  Wet_Air
                     icetray.OMKey(79, 32),  # TP9P3835  Drift
                     icetray.OMKey(79, 41),  # UP9P3832  Flash_Flood
                     icetray.OMKey(79, 53),  # UP9P3834  Temperature_Inversion
                     icetray.OMKey(79, 60),  # TP9P3905  General_Wetness
                     icetray.OMKey(14,  3),  # UP9P3902  Gods_Tears
                     icetray.OMKey(14,  8),  # TP9P3853  Atlas_Bear
                     icetray.OMKey(14, 14),  # TP9P3859  Anancus
                     icetray.OMKey(14, 21),  # UP9P3914  High_Fog
                     icetray.OMKey(14, 28),  # TP9P3909  Gullywasher
                     icetray.OMKey(14, 41),  # UP9P3900  Ghostly_Fog
                     icetray.OMKey(14, 51),  # UP9P3894  Spitting
                     icetray.OMKey(14, 58)]) # TP9P3849  Cretan_Giant_Owl

    # LED colors by flasher mask index for cDOMs
    # (from http://wiki.icecube.wisc.edu/index.php/CDOM_Info)
    colorDOMledColors = [I3CLSimFlasherPulse.FlasherPulseType.LED505nm,  # LED1  (0)      (narrow beam)
                         I3CLSimFlasherPulse.FlasherPulseType.LED450nm,  # LED2  (1)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED505nm,  # LED3  (2)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED450nm,  # LED4  (3)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED505nm,  # LED5  (4)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED450nm,  # LED6  (5)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED340nm,  # LED7  (6)      ( wide beam )
                         I3CLSimFlasherPulse.FlasherPulseType.LED370nm,  # LED8  (7)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED340nm,  # LED9  (8)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED370nm,  # LED10 (9)      (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED340nm,  # LED11 (10)     (     "     )
                         I3CLSimFlasherPulse.FlasherPulseType.LED370nm]  # LED12 (11)     (     "     )

    # LED angular emission profile by color
    #
    # (gaussian sigma in polar ([0]) and azimuthal ([1]) direction)
    # Measured in air and converted to ice using a MC simulation for 405nm.
    # Other wavelengths were measured in air and re-scaled according to the
    # 405nm scaling w/o a dedicated MC simulation.
    #
    # (from http://wiki.icecube.wisc.edu/index.php/CDOM_Info and
    # http://wiki.icecube.wisc.edu/index.php/LED_angular_emission_profile)
    LEDangularEmissionProfile = {(I3CLSimFlasherPulse.FlasherPulseType.LED405nm, True)  : [ 9.7*I3Units.deg,  9.8*I3Units.deg], # tilted
                                 (I3CLSimFlasherPulse.FlasherPulseType.LED405nm, False) : [ 9.2*I3Units.deg, 10.1*I3Units.deg], # horizontal
                                 # cDOM LEDs are all horizontal
                                 (I3CLSimFlasherPulse.FlasherPulseType.LED340nm, False) : [36.1*I3Units.deg, 39.6*I3Units.deg], # the wiki page says 42.9 instead of 39.6, but this seems to be a copy&paste error
                                 (I3CLSimFlasherPulse.FlasherPulseType.LED370nm, False) : [39.1*I3Units.deg, 42.9*I3Units.deg],
                                 (I3CLSimFlasherPulse.FlasherPulseType.LED450nm, False) : [ 4.8*I3Units.deg,  5.3*I3Units.deg],
                                 (I3CLSimFlasherPulse.FlasherPulseType.LED505nm, False) : [ 4.5*I3Units.deg,  4.9*I3Units.deg]}

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("FlasherInfoVectName",
                          "Name of the I3FlasherInfoVect object to read",
                          "I3FlasherInfo")
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3CLSimFlasherPulseSeries to write",
                          "I3CLSimFlasherPulseSeries")
        self.AddParameter("NumberOfPhotonsAtMaxBrightness",
                          "Number of photons generated by a flasher LED at maximum brightness",
                          1.17e10) # this number has been determined by Dima from SPICE-Lea ice model fits
        self.AddParameter("FlasherOMKeyVectName",
                          "Name of a vector of OMKeys with all active flasher OMs",
                          None)
        self.AddOutBox("OutBox")
        self.currentGeometry = None
        
    def Configure(self):
        self.flasherInfoVectName = self.GetParameter("FlasherInfoVectName")
        self.flasherPulseSeriesName = self.GetParameter("FlasherPulseSeriesName")
        self.numberOfPhotonsAtMaxBrightness = self.GetParameter("NumberOfPhotonsAtMaxBrightness")
        self.flasherOMKeyVectName = self.GetParameter("FlasherOMKeyVectName")

    def Geometry(self, frame):
        self.currentGeometry = frame["I3Geometry"]
        self.PushFrame(frame)
    
    def GetNumPhotons(self, brightness, width):
        # taken from Chris Wendt's wiki page http://wiki.icecube.wisc.edu/index.php/LED_light_output
        # (this is for a single LED)
        return self.numberOfPhotonsAtMaxBrightness*(0.0006753 + 0.00005593 * float(brightness)) * (float(width) + 13.9 - (57.5/(1. + float(brightness) / 34.4)))
    
    def DAQ(self, frame):
        if self.currentGeometry is None:
            raise RuntimeError("found a DAQ frame without any previous Geometry frame")
        
        outputSeries = I3CLSimFlasherPulseSeries()
        if self.flasherInfoVectName not in frame:
            # put an empty pulseSeries in the frame
            frame[self.flasherPulseSeriesName] = outputSeries
            self.PushFrame(frame)
            return
        
        flasherOMKeyVect = None
        if self.flasherOMKeyVectName is not None and self.flasherOMKeyVectName!="":
            flasherOMKeyVect = dataclasses.I3VectorOMKey()
        
        flasherInfoVect = frame[self.flasherInfoVectName] 
        for flasherInfo in flasherInfoVect:
            if flasherOMKeyVect is not None:
                flasherOMKeyVect.append(flasherInfo.flashing_om)
            
            isColorDOM = flasherInfo.flashing_om in FlasherInfoVectToFlasherPulseSeriesConverter.colorDOMs
            
            omGeo = self.currentGeometry.omgeo[flasherInfo.flashing_om]
            domPos = omGeo.position
            
            numPhotons = self.GetNumPhotons(flasherInfo.led_brightness, flasherInfo.width)
            
            # loop over all twelve flashers
            for i in range(12):
                # skip LED if inactive (check in flasher mask)
                if flasherInfo.mask & (1 << i) == 0: continue
                
                # is the LED tilted?
                if isColorDOM:
                    tiltedFlasher = False # cDOM flashers are all horizontal
                else:
                    tiltedFlasher = i < 6 # 0-5 are tilted, 6-11 are horizontal
                    
                flasherPosIndex = i % 6
                flasherRadiusInDOM = 11.9*I3Units.cm
                flasherZPosInDOM = 8.0*I3Units.cm

                if "I3Orientation" in dir(dataclasses):
                    # new-style I3OMGeo (with full DOM orientation information):
                    # do calculations in DOM coordinate frame first (i.e. DOM axis points towards +z)

                    # the rotation of flashers is counter-clockwise when viewed from above
                    # (in DOM frame where the DOM points upwards).
                    # Flasher LED1 and LED7 (both flasherPosIndex==0) are assumed 
                    # to point towards the DOM's x-axis.
                    flasherAziRotation = 60.*I3Units.deg * float(flasherPosIndex)
                    flasherUpwardsTiltAngle = 0.
                    if tiltedFlasher: flasherUpwardsTiltAngle=48.*I3Units.deg # (or 42deg down from the vertical)

                    flasherDir = dataclasses.I3Direction()
                    flasherDir.set_theta_phi(90.*I3Units.deg + flasherUpwardsTiltAngle, flasherAziRotation)
                    # now apply the DOM orientation (this will flip over the DOM)
                    flasherDir = omGeo.orientation.rotate_out(flasherDir)
                    #flasherDir = omGeo.orientation.rotate(flasherDir)
                    # for standard DOMs (i.e. pointing down), flashers are now 
                    # ordered clockwise when viewed from above
                    
                    flasherPos = dataclasses.I3Direction()
                    # flasher position on flasher board
                    flasherPos.set_theta_phi(90.*I3Units.deg, flasherAziRotation) # no upwards tilt for position
                    # flasher position in DOM coordinates
                    flasherPos = dataclasses.I3Position(flasherPos.x*flasherRadiusInDOM, flasherPos.y*flasherRadiusInDOM, flasherZPosInDOM)
                    # apply the DOM orientation (this will flip over the DOM)
                    #flasherPos = omGeo.orientation.rotate(flasherPos)
                    flasherPos = omGeo.orientation.rotate_out(flasherPos)
                    # shift to the actual position in global coordinates
                    flasherPos = dataclasses.I3Position(domPos.x + flasherPos.x, domPos.y + flasherPos.y, domPos.z - flasherPos.z)
                else:
                    # old-style I3OMGeo (only up/down enum and aziangle for DOM orientation):
                    # do calculations in global coordinate system (DOM axis towards -z)

                    # the rotation of flashers is clockwise when viewed from above
                    # flasher LED1 and LED7 (both flasherPosIndex==0) are assumed 
                    # to point towards the DOM's x-axis.
                    flasherAziRotation = -60.*I3Units.deg * float(flasherPosIndex)
                    flasherUpwardsTiltAngle = 0.
                    if tiltedFlasher: flasherUpwardsTiltAngle=48.*I3Units.deg # (or 42deg down from the vertical)

                    flasherDir = dataclasses.I3Direction()
                    flasherDir.set_theta_phi(90.*I3Units.deg - flasherUpwardsTiltAngle, flasherAziRotation - omGeo.aziangle)

                    flasherPos = dataclasses.I3Direction()
                    # flasher position on flasher board
                    flasherPos.set_theta_phi(90.*I3Units.deg, flasherAziRotation - omGeo.aziangle) # no upwards tilt for position
                    # flasher position in global coordinates
                    flasherPos = dataclasses.I3Position(domPos.x + flasherPos.x*flasherRadiusInDOM, domPos.y + flasherPos.y*flasherRadiusInDOM, domPos.z + flasherZPosInDOM)

                newPulse = I3CLSimFlasherPulse()

                if isColorDOM:
                    newPulse.type = FlasherInfoVectToFlasherPulseSeriesConverter.colorDOMledColors[i]
                else:
                    newPulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
                newPulse.pos = flasherPos
                newPulse.dir = flasherDir
                newPulse.time = flasherInfo.flash_time
                # FWHM according to http://wiki.icecube.wisc.edu/index.php/LED_output_time_profile :
                newPulse.pulseWidth = (float(flasherInfo.width)/2.) * I3Units.ns
                newPulse.numberOfPhotonsNoBias = numPhotons

                newPulse.angularEmissionSigmaPolar = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(newPulse.type, tiltedFlasher)][0]
                newPulse.angularEmissionSigmaAzimuthal = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(newPulse.type, tiltedFlasher)][1]
                
                outputSeries.append(newPulse)

        frame[self.flasherPulseSeriesName] = outputSeries
        
        if flasherOMKeyVect is not None:
            frame[self.flasherOMKeyVectName] = flasherOMKeyVect
        
        self.PushFrame(frame)
    
