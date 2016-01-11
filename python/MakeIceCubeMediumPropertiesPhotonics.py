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
# $Id: MakeIceCubeMediumPropertiesPhotonics.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file MakeIceCubeMediumPropertiesPhotonics.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from __future__ import print_function

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimMediumProperties
from icecube.clsim import I3CLSimRandomValueMixed
from icecube.clsim import I3CLSimRandomValueHenyeyGreenstein
from icecube.clsim import I3CLSimRandomValueSimplifiedLiu
from icecube.clsim import I3CLSimFunctionFromTable
from icecube.clsim import I3CLSimFunctionRefIndexIceCube
from icecube.clsim import I3CLSimScalarFieldConstant
from icecube.clsim import I3CLSimVectorTransformConstant

from I3Tray import I3Units

import numpy, math
from os.path import expandvars


def MakeIceCubeMediumPropertiesPhotonics(tableFile,
                                         detectorCenterDepth = 1948.07*I3Units.m):
    ### read ice information from a photonics-compatible table

    file = open(tableFile, 'r')
    rawtable = file.readlines()
    file.close()
    del file

    # get rid of comments and split lines at whitespace
    parsedtable = [line.split() for line in rawtable if line.lstrip()[0] != "#"]

    nlayer_lines = [line for line in parsedtable if line[0].upper()=="NLAYER"]
    nwvl_lines =   [line for line in parsedtable if line[0].upper()=="NWVL"]

    if len(nlayer_lines) == 0:
        raise RuntimeError("There is no \"NLAYER\" entry in your ice table!")
    if len(nlayer_lines) > 1:
        raise RuntimeError("There is more than one \"NLAYER\" entry in your ice table!")
    if len(nwvl_lines) == 0:
        raise RuntimeError("There is no \"NWVL\" entry in your ice table!")
    if len(nwvl_lines) > 1:
        raise RuntimeError("There is more than one \"NWVL\" entry in your ice table!")

    nLayers = int(nlayer_lines[0][1])
    nWavelengths = int(nwvl_lines[0][1])
    startWavelength = float(nwvl_lines[0][2])*I3Units.nanometer
    stepWavelength = float(nwvl_lines[0][3])*I3Units.nanometer

    startWavelength += stepWavelength/2.

    # make a list of wavelengths
    wavelengths = numpy.linspace(start=startWavelength, stop=startWavelength+float(nWavelengths)*stepWavelength, num=nWavelengths, endpoint=False)

    print("The table file {0} has {1} layers and {2} wavelengths, starting from {3}ns in {4}ns steps".format(tableFile, nLayers, nWavelengths, startWavelength/I3Units.nanometer, stepWavelength/I3Units.nanometer))


    # replace parsedtable with a version without NLAYER and NWVL entries
    parsedtable = [line for line in parsedtable if ((line[0].upper()!="NLAYER") and (line[0].upper()!="NWVL"))]

    if len(parsedtable) != nLayers*6:
        raise RuntimeError("Expected {0}*6={1} lines [not counting NLAYER and NWVL] in the icetable file, found {2}".format(nLayers, nLayers*6, len(parsedtable)))

    if parsedtable[0][0].upper() != "LAYER":
        raise RuntimeError("Layer definitions should start with the LAYER keyword (reading {0})".format(tableFile))

    # parse the lines
    layers = []

    layerdict = dict()
    for line in parsedtable:
        keyword = line[0].upper()

        if keyword == "LAYER":
            if len(layerdict) != 0:
                layers.append(layerdict)

            #re-set layerdict
            layerdict = dict()
        else:
            # some other keyword
            if keyword in layerdict:
                raise RuntimeError("Keyword {0} is used twice for one layer (reading {1})".format(keyword, tableFile))

        # optional units
        unit = 1.                                       # no units by default
        if keyword is "LAYER": unit = I3Units.meter     # coordinates
        if keyword is "ABS": unit = 1./I3Units.meter    # absorption coefficient
        if keyword is "SCAT": unit = 1./I3Units.meter   # scattering coefficient

        layerdict[keyword] = [float(number)*unit for number in line[1:]]

    # last layer
    if len(layerdict) != 0:
        layers.append(layerdict)

    if len(layers) < 1:
        raise RuntimeError("At least one layers is requried (reading {0})".format(tableFile))

    layerHeight = abs(layers[0]['LAYER'][1] - layers[0]['LAYER'][0])
    print("layer height is {0}m".format(layerHeight/I3Units.m))

    # sort layers into dict by bottom z coordinate and check some assumptions
    layersByZ = dict()
    for layer in layers:
        layerBottom = layer['LAYER'][0]
        layerTop = layer['LAYER'][1]

        if layerBottom > layerTop:
            print("a layer is upside down. compensating. (reading {0})".format(tableFile))
            dummy = layerBottom
            layerBottom = layerTop
            layerTop = dummy

        # check layer height
        if abs((layerTop-layerBottom) - layerHeight) > 0.0001:
            raise RuntimeError("Differing layer heights while reading {0}. Expected {1}m, got {2}m.".format(tableFile, layerHeight, layerTop-layerBottom))

        layersByZ[layerBottom] = layer


    layers = []
    endZ = None
    for dummy, layer in sorted(layersByZ.items()):
        startZ = layer['LAYER'][0]
        if (endZ is not None) and (abs(endZ-startZ) > 0.0001):
            raise RuntimeError("Your layers have holes. previous layer ends at {0}m, next one starts at {1}m".format(endZ/I3Units.m, startZ/I3Units.m))
        endZ = layer['LAYER'][1]
        layers.append(layer)

    # check even more assumptions
    meanCos=None
    for layer in layers:
        if meanCos is None: meanCos = layer['COS'][0]

        for cosVal in layer['COS']:
            if abs(cosVal-meanCos) > 0.0001:
                raise RuntimeError("only a constant mean cosine is supported by clsim (expected {0}, got {1})".format(meanCos, cosVal))

        if len(layer['COS']) != len(wavelengths):
            raise RuntimeError("Expected {0} COS values, got {1}".format(len(wavelengths), len(layer['COS'])))
        if len(layer['ABS']) != len(wavelengths):
            raise RuntimeError("Expected {0} ABS values, got {1}".format(len(wavelengths), len(layer['ABS'])))
        if len(layer['SCAT']) != len(wavelengths):
            raise RuntimeError("Expected {0} SCAT values, got {1}".format(len(wavelengths), len(layer['SCAT'])))
        if len(layer['N_GROUP']) != len(wavelengths):
            raise RuntimeError("Expected {0} N_GROUP values, got {1}".format(len(wavelengths), len(layer['N_GROUP'])))
        if len(layer['N_PHASE']) != len(wavelengths):
            raise RuntimeError("Expected {0} N_PHASE values, got {1}".format(len(wavelengths), len(layer['N_PHASE'])))

        for i in range(len(wavelengths)):
            if abs(layer['N_GROUP'][i]-layers[0]['N_GROUP'][i]) > 0.0001:
                raise RuntimeError("N_GROUP may not be different for different layers in this version of clsim!")

            if abs(layer['N_PHASE'][i]-layers[0]['N_PHASE'][i]) > 0.0001:
                raise RuntimeError("N_PHASE may not be different for different layers in this version of clsim!")

        #print "start", layer['LAYER'][0], "end", layer['LAYER'][1]


    ##### start making the medium property object

    m = I3CLSimMediumProperties(mediumDensity=0.9216*I3Units.g/I3Units.cm3,
                                layersNum=len(layers),
                                layersZStart=layers[0]['LAYER'][0],
                                layersHeight=layerHeight,
                                rockZCoordinate=-870.*I3Units.m,
                                # TODO: inbetween: from 1740 upwards: less dense ice (factor 0.825)
                                airZCoordinate=1940.*I3Units.m)


    iceCubeScatModel = I3CLSimRandomValueHenyeyGreenstein(meanCosine=meanCos)
    m.SetScatteringCosAngleDistribution(iceCubeScatModel)

    # no ice/water anisotropy. all of these three are no-ops
    m.SetDirectionalAbsorptionLengthCorrection(I3CLSimScalarFieldConstant(1.))
    m.SetPreScatterDirectionTransform(I3CLSimVectorTransformConstant())
    m.SetPostScatterDirectionTransform(I3CLSimVectorTransformConstant())

    # no ice tilt
    m.SetIceTiltZShift(I3CLSimScalarFieldConstant(0.))

    phaseRefIndex = I3CLSimFunctionFromTable(startWavelength, stepWavelength, layers[0]['N_PHASE'])
    groupRefIndex = I3CLSimFunctionFromTable(startWavelength, stepWavelength, layers[0]['N_GROUP'])
    #phaseRefIndex = I3CLSimFunctionRefIndexIceCube(mode="phase")
    #groupRefIndex = I3CLSimFunctionRefIndexIceCube(mode="group")

    for i in range(len(layers)):
        #print "layer {0}: depth at bottom is {1} (z_bottom={2}), b_400={3}".format(i, depthAtBottomOfLayer[i], layerZStart[i], b_400[i])

        m.SetPhaseRefractiveIndex(i, phaseRefIndex)
        m.SetGroupRefractiveIndexOverride(i, groupRefIndex)

        absLenTable = [1./absCoeff for absCoeff in layers[i]['ABS']]
        absLen = I3CLSimFunctionFromTable(startWavelength, stepWavelength, absLenTable, storeDataAsHalfPrecision=True)
        m.SetAbsorptionLength(i, absLen)

        scatLenTable = [(1./scatCoeff)*(1.-meanCos) for scatCoeff in layers[i]['SCAT']]
        scatLen = I3CLSimFunctionFromTable(startWavelength, stepWavelength, scatLenTable, storeDataAsHalfPrecision=True)
        m.SetScatteringLength(i, scatLen)

    return m

