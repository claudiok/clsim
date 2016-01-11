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
# $Id: GetIceTiltZShift.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetIceTiltZShift.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from os.path import expandvars

import numpy
from icecube.clsim import I3CLSimScalarFieldIceTiltZShift
from icecube.icetray import I3Units

try:
    from icecube.dataclasses import I3Matrix
except ImportError:
    # I3Matrix is in dataclasses trunk, but not in
    # an offline-software release yet.. :-(
    from icecube.clsim import I3Matrix

def GetIceTiltZShift(
    tiltDirAzimuth = 225.*I3Units.deg,
    tiltDirectory = expandvars("$I3_BUILD/clsim/resources/ice/TILT_data/"),
    detectorCenterDepth = 1948.07*I3Units.m
    ):

    distance_from_origin_in_tilt_dir = numpy.loadtxt(tiltDirectory + "/tilt.par", unpack=True)[1]*I3Units.m

    tilt_dat = numpy.loadtxt(tiltDirectory + "/tilt.dat", unpack=True)
    zcoords = (detectorCenterDepth-tilt_dat[0])[::-1]

    zshift = []
    for i in range(len(distance_from_origin_in_tilt_dir)):
        zshift.append(tilt_dat[i+1][::-1])
    zshift = numpy.array(zshift)

    return I3CLSimScalarFieldIceTiltZShift(
        distancesFromOriginAlongTilt = distance_from_origin_in_tilt_dir,
        zCoordinates = zcoords,
        zCorrections = I3Matrix(zshift),
        directionOfTiltAzimuth = tiltDirAzimuth)


