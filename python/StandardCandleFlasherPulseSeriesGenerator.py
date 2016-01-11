#
# Copyright (c) 2012
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
# $Id: StandardCandleFlasherPulseSeriesGenerator.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file StandardCandleFlasherPulseSeriesGenerator.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#


from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries

from I3Tray import I3Units

class StandardCandleFlasherPulseSeriesGenerator(icetray.I3ConditionalModule):
    """
    Generates I3CLSimFlasherPulse objects for IceCube Standard Candle (I&II)
    simulation.
    
    Values are taken from internal report icecube/200704001:
    http://internal.icecube.wisc.edu/reports/data/icecube/2007/04/001/icecube_200704001_v1.pdf
    """
    
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("FlasherPulseSeriesName",
                          "Name of the I3CLSimFlasherPulseSeries to write",
                          "I3CLSimFlasherPulseSeries")
        self.AddParameter("PhotonsPerPulse",
                          "Photons to emit per pulse (see http://wiki.icecube.wisc.edu/index.php/Standard_Candle\n" +
                          "for some more information on what number to use)",
                          2.5e13)
        self.AddParameter("FlashTime",
                          "Time (within each event) at which to flash",
                          0.*I3Units.ns)
        self.AddParameter("CandleNumber",
                          "Choose which source to simulate: 1 for SC1 or 2 for SC2",
                          1)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.flasherPulseSeriesName = self.GetParameter("FlasherPulseSeriesName")
        self.photonsPerPulse = self.GetParameter("PhotonsPerPulse")
        self.flashTime = self.GetParameter("FlashTime")
        self.candleNumber = self.GetParameter("CandleNumber")
        
        if self.candleNumber not in [1,2]:
            raise RuntimeError("You have to select either standard candle 1 or 2. You chose %u." % self.candleNumber)

    def DAQ(self, frame):
        outputSeries = I3CLSimFlasherPulseSeries()

        numPhotons = self.photonsPerPulse

        newPulse = I3CLSimFlasherPulse()
        if self.candleNumber==1:
            newPulse.type = I3CLSimFlasherPulse.FlasherPulseType.SC1
            # from http://wiki.icecube.wisc.edu/index.php/File:Sc_geometry.jpg
            newPulse.pos = dataclasses.I3Position(544.07*I3Units.m, 55.89*I3Units.m, 136.86*I3Units.m)
            newPulse.dir = dataclasses.I3Direction(0.,0., 1.) # facing up
        elif self.candleNumber==2:
            newPulse.type = I3CLSimFlasherPulse.FlasherPulseType.SC2
            # from http://wiki.icecube.wisc.edu/index.php/Standard_Candle#Standard_Candle_Mark_II
            newPulse.pos = dataclasses.I3Position(11.87*I3Units.m, 179.19*I3Units.m, -205.64*I3Units.m)
            newPulse.dir = dataclasses.I3Direction(0.,0.,-1.) # facing down
        else:
            raise RuntimeError("invalid candle number. logic error.")

        newPulse.time = self.flashTime
        newPulse.numberOfPhotonsNoBias = self.photonsPerPulse

        # from icecube/200704001 section 2:
        newPulse.pulseWidth = 4. * I3Units.ns

        # these two have different meanings than for flashers:
        # polar is the angle w.r.t. the candle axis,
        # azimuthal is angle along the circle (and should always be 360deg)
        newPulse.angularEmissionSigmaPolar = 41.13*I3Units.deg
        newPulse.angularEmissionSigmaAzimuthal = 360.*I3Units.deg

        # insert a single pulse
        outputSeries.append(newPulse)

        frame[self.flasherPulseSeriesName] = outputSeries

        self.PushFrame(frame)
