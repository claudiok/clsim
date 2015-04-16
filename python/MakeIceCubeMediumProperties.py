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
# $Id$
# 
# @file MakeIceCubeMediumProperties.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

from __future__ import print_function

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimMediumProperties
from icecube.clsim import I3CLSimRandomValueMixed
from icecube.clsim import I3CLSimRandomValueHenyeyGreenstein
from icecube.clsim import I3CLSimRandomValueSimplifiedLiu
from icecube.clsim import I3CLSimFunctionRefIndexIceCube
from icecube.clsim import I3CLSimFunctionAbsLenIceCube
from icecube.clsim import I3CLSimFunctionScatLenIceCube
from icecube.clsim import I3CLSimScalarFieldConstant
from icecube.clsim import I3CLSimVectorTransformConstant

from . import util

from I3Tray import I3Units

import numpy, math
import os
from os.path import expandvars


def MakeIceCubeMediumProperties(detectorCenterDepth = 1948.07*I3Units.m,
                                iceDataDirectory=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                                useTiltIfAvailable=True,
                                returnParameters=False):
    ### read ice information from PPC-compatible tables
    
    # do we have tilt descripton files?
    useTilt=False
    if useTiltIfAvailable:
        hasTiltPar = os.path.isfile(iceDataDirectory+"/tilt.par")
        hasTiltDat = os.path.isfile(iceDataDirectory+"/tilt.dat")
        if hasTiltPar and not hasTiltDat:
            raise RuntimeError("ice model directory has tilt.par but tilt.dat is missing!")
        elif hasTiltDat and not hasTiltPar:
            raise RuntimeError("ice model directory has tilt.dat but tilt.par is missing!")
        elif hasTiltDat and hasTiltPar:
            useTilt=True

    icemodel_dat = numpy.loadtxt(iceDataDirectory+"/icemodel.dat", unpack=True)
    icemodel_par = numpy.loadtxt(iceDataDirectory+"/icemodel.par")
    icemodel_cfg = numpy.loadtxt(iceDataDirectory+"/cfg.txt")
    
    # is this Dima's new 4-parameter file or his old 6-parameter file?
    if len(icemodel_par)==6:
        alpha = icemodel_par[0][0]
        kappa = icemodel_par[1][0]
        A     = icemodel_par[2][0]
        B     = icemodel_par[3][0]
        D     = icemodel_par[4][0]
        E     = icemodel_par[5][0]
    elif len(icemodel_par)==4:
        # the 4-parameter files appeared up in ppc around March 2012
        alpha = icemodel_par[0][0]
        kappa = icemodel_par[1][0]
        A     = icemodel_par[2][0]
        B     = icemodel_par[3][0]
        
        # this is what ppc does to fill the remaining two parameters:
        wv0 = 400.
        D     = wv0**kappa
        E     = 0.
    else:
        raise RuntimeError(iceDataDirectory+"/icemodel.par is not a valid Dima-icemodel file. (needs either 4 or 6 entries, this one has %u entries)" % len(icemodel_par))
    
    if len(icemodel_cfg) < 4:
        raise RuntimeError(iceDataDirectory+"/cfg.txt does not have enough configuration lines. It needs at least 4.")

    oversizeScaling       = icemodel_cfg[0] # currently ignored
    efficiencyCorrection  = icemodel_cfg[1] # currently ignored
    liuScatteringFraction = icemodel_cfg[2]
    meanCosineTheta       = icemodel_cfg[3]

    hasAnisotropy = False
    if len(icemodel_cfg) > 4 and len(icemodel_cfg) < 7:
        raise RuntimeError(iceDataDirectory+"/cfg.txt has more than 4 lines (this means you probably get ice anisotropy), but it needs at least 7 lines in this case.")
    elif len(icemodel_cfg) > 4:
        hasAnisotropy = True
        anisotropyDirAzimuth  = icemodel_cfg[4]*I3Units.deg # direction of ice tilt (perp. to flow)
        magnitudeAlongDir     = icemodel_cfg[5]             # magnitude of ice anisotropy along tilt
        magnitudePerpToDir    = icemodel_cfg[6]             # magnitude of ice anisotropy along flow


    if liuScatteringFraction<0. or liuScatteringFraction>1.:
        raise RuntimeError("Invalid Liu(SAM) scattering fraction configured in cfg.txt: value=%g" % liuScatteringFraction)
    if meanCosineTheta<-1. or meanCosineTheta>1.:
        raise RuntimeError("Invalid <cos(theta)> configured in cfg.txt: value=%g" % meanCosineTheta)
    
    depth = icemodel_dat[0]*I3Units.m
    b_e400 = icemodel_dat[1]
    a_dust400 = icemodel_dat[2]
    delta_tau = icemodel_dat[3]
    
    # check delta_tau values against formula
    # According to the IceCube paper, these values should be calculated like this.
    # Values in tables seem to have an offset of 6.53m
    def temp(depth):
        return 221.5 - 0.00045319*(depth/I3Units.m) + 5.822e-6 * (depth/I3Units.m)**2.
    def deltaTau(depth):
        return temp(depth+6.53*I3Units.m) - temp(1730.*I3Units.m)
    
    maxRelativeError=0.
    for thisDepth, thisDeltaTau in numpy.array([depth, delta_tau]).transpose():
        relativeError = abs(thisDeltaTau-deltaTau(thisDepth))/thisDeltaTau
        if relativeError > maxRelativeError: maxRelativeError=relativeError
    if maxRelativeError > 0.01:
        print("The ice table's delta_tau values do not seem to correpsond to the equation. Loading table anyway.")
    
    
    # some sanity checks on the layer files
    if len(depth)<2: raise RuntimeError("There is only a single layer in your layer definition file")
    
    layerHeight = depth[1]-depth[0]
    if layerHeight<=0.: raise RuntimeError("ice layer depths are not in increasing order")
    
    for i in range(len(depth)-1):
        if abs((depth[i+1]-depth[i]) - layerHeight) > 1e-5:
            raise RuntimeError("ice layers are not spaced evenly")
    
    # change the order (top-to-bottom -> bottom-to-top)
    depth = depth[::-1]
    b_e400 = b_e400[::-1]               # effective scattering length
    a_dust400 = a_dust400[::-1]
    delta_tau = delta_tau[::-1]

    b_400 = b_e400/(1.-meanCosineTheta) # scattering length (used in the simulation)
    
    # Be compatible with PPC, which assumes the specified depths
    # are in the middle of the layer. We need the depth at the
    # top of the layer here, so correct for that:
    depth = depth-layerHeight/2.
    
    # layerZ is in z-coordinates, from bottom to top (ascending z)
    depthAtBottomOfLayer = depth + layerHeight
    layerZStart = detectorCenterDepth - depthAtBottomOfLayer
    layerZEnd = detectorCenterDepth - depth
    
    ##### start making the medium property object
    
    m = I3CLSimMediumProperties(mediumDensity=0.9216*I3Units.g/I3Units.cm3,
                                layersNum=len(layerZStart),
                                layersZStart=layerZStart[0],
                                layersHeight=layerHeight,
                                rockZCoordinate=-870.*I3Units.m,
                                # TODO: inbetween: from 1740 upwards: less dense ice (factor 0.825)
                                airZCoordinate=1940.*I3Units.m)
    
    # None of the IceCube wlen-dependent functions have a fixed minimum
    # or maximum wavelength value set. We need to set some sensible range here.
    # This seems to be the definition range of the DOM wavelength acceptance:
    m.ForcedMinWlen = 265.*I3Units.nanometer
    m.ForcedMaxWlen = 675.*I3Units.nanometer
    
    iceCubeScatModel = I3CLSimRandomValueMixed(
        firstDistribution=I3CLSimRandomValueSimplifiedLiu(meanCosine=meanCosineTheta),
        secondDistribution=I3CLSimRandomValueHenyeyGreenstein(meanCosine=meanCosineTheta),
        fractionOfFirstDistribution=liuScatteringFraction)
    m.SetScatteringCosAngleDistribution(iceCubeScatModel)
    
    if not hasAnisotropy:
        # no ice/water anisotropy. all of these three are no-ops
        m.SetDirectionalAbsorptionLengthCorrection(I3CLSimScalarFieldConstant(1.))
        m.SetPreScatterDirectionTransform(I3CLSimVectorTransformConstant())
        m.SetPostScatterDirectionTransform(I3CLSimVectorTransformConstant())
    else:
        # print("Anisotropy! Whooo!", anisotropyDirAzimuth/I3Units.deg, magnitudeAlongDir, magnitudePerpToDir)

        absLenScaling, preScatterTransform, postScatterTransform = \
            util.GetSpiceLeaAnisotropyTransforms(
                anisotropyDirAzimuth,
                magnitudeAlongDir,
                magnitudePerpToDir
                )

        m.SetDirectionalAbsorptionLengthCorrection(absLenScaling)
        m.SetPreScatterDirectionTransform(preScatterTransform)
        m.SetPostScatterDirectionTransform(postScatterTransform)

    if useTilt:
        # print("Tilt! Wheee!")

        m.SetIceTiltZShift(
            util.GetIceTiltZShift(
                tiltDirectory = iceDataDirectory,
                detectorCenterDepth = detectorCenterDepth,
                )
            )
    else:
        # no ice tilt
        m.SetIceTiltZShift(I3CLSimScalarFieldConstant(0.))

    phaseRefIndex = I3CLSimFunctionRefIndexIceCube(mode="phase")
    groupRefIndex = I3CLSimFunctionRefIndexIceCube(mode="group")
    for i in range(len(layerZStart)):
        #print "layer {0}: depth at bottom is {1} (z_bottom={2}), b_400={3}".format(i, depthAtBottomOfLayer[i], layerZStart[i], b_400[i])
        
        m.SetPhaseRefractiveIndex(i, phaseRefIndex)

        # the IceCube group refractive index parameterization is not exactly 
        # what you would expect from the phase refractive index.
        # To use calculated values instead of the parameterization,
        # just comment this line:
        m.SetGroupRefractiveIndexOverride(i, groupRefIndex)
        
        absLen = I3CLSimFunctionAbsLenIceCube(kappa=kappa, A=A, B=B, D=D, E=E,
                                                        aDust400=a_dust400[i],
                                                        deltaTau=delta_tau[i])
        m.SetAbsorptionLength(i, absLen)

        scatLen = I3CLSimFunctionScatLenIceCube(alpha=alpha,
                                                          b400=b_400[i])
        m.SetScatteringLength(i, scatLen)

    if not returnParameters:
        return m
    else:
        parameters = dict()
        if hasAnisotropy:
            parameters["anisotropyDirAzimuth"]=anisotropyDirAzimuth
            parameters["anisotropyMagnitudeAlongDir"]=magnitudeAlongDir
            parameters["anisotropyMagnitudePerpToDir"]=magnitudePerpToDir
        else:
            parameters["anisotropyDirAzimuth"]=float('NaN')
            parameters["anisotropyMagnitudeAlongDir"]=float('NaN')
            parameters["anisotropyMagnitudePerpToDir"]=float('NaN')
        return (m, parameters)
