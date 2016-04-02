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
# $Id: GetSpiceLeaAnisotropyTransforms.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetSpiceLeaAnisotropyTransforms.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

import numpy
from icecube.clsim import I3CLSimScalarFieldAnisotropyAbsLenScaling
from icecube.clsim import I3CLSimVectorTransformMatrix
from icecube.icetray import I3Units

try:
    from icecube.dataclasses import I3Matrix
except ImportError:
    # I3Matrix is in dataclasses trunk, but not in
    # an offline-software release yet.. :-(
    from icecube.clsim import I3Matrix

def GetSpiceLeaAnisotropyTransforms(
    anisotropyDirAzimuth=216.*I3Units.deg, # direction of ice tilt (perp. to flow)
    magnitudeAlongDir=0.04,                # magnitude of ice anisotropy along tilt
    magnitudePerpToDir=-0.08):             # magnitude of ice anisotropy along flow
    """
    Returns the direction-dependent absorption length
    scaling function and a pre- and post-scattering
    direction transformation used in the ice anisotropy
    description of Spice-Lea.
    """
    
    absLenScaling = I3CLSimScalarFieldAnisotropyAbsLenScaling(
        anisotropyDirAzimuth=anisotropyDirAzimuth,
        magnitudeAlongDir=magnitudeAlongDir,
        magnitudePerpToDir=magnitudePerpToDir)

    # The matrix A in Dima's description is the distortion
    # in the frame of (n1,n2,n3), where n1 points in the
    # direction of the anisotropy and n2 points into the
    # direction perpendicular to it.
    k1 = numpy.exp(magnitudeAlongDir)  # coefficient of anisotropy parallel to "tilt" direction
    k2 = numpy.exp(magnitudePerpToDir) # coefficient of anisotropy perpendicular to "tilt" direction
    kz = 1./(k1*k2)                    # a normalizing factor for the z direction
    A = numpy.array(
        [[k1, 0., 0.],
         [0., k2, 0.],
         [0., 0., kz]]
        )

    # The orthogonal matrix T transforms a vector into
    # the (n1,n2,n3) coordinate system. (It is a simple
    # rotation around the z-axis.)
    sa = numpy.sin(anisotropyDirAzimuth)
    ca = numpy.cos(anisotropyDirAzimuth)
    T = numpy.array(
        [[ ca, sa, 0.],
         [-sa, ca, 0.],
         [ 0., 0., 1.]]
        )

    # Now combine the transformations: rotate into the
    # new coordinate system first, apply A and rotate back.
    # This is applied before the photon is rotated into
    # the new direction.
    # (read this from right to left)
    Cpre = numpy.dot(numpy.dot(T.T, A), T)

    preScatterTransform = I3CLSimVectorTransformMatrix(
        I3Matrix(Cpre),
        renormalize=True) # after applying Cpost, make it a unit vector again

    # Now combine the transformations: rotate into the
    # new coordinate system first, apply A^-1 and rotate back.
    # This is applied after the photon has been rotated into
    # its new direction.
    # (read this from right to left)
    Cpost = numpy.dot(numpy.dot(T.T, numpy.linalg.inv(A)), T)

    postScatterTransform = I3CLSimVectorTransformMatrix(
        I3Matrix(Cpost),
        renormalize=True) # after applying Cpost, make it a unit vector again

    return (absLenScaling, preScatterTransform, postScatterTransform)


