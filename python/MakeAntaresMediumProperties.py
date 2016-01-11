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
# $Id: MakeAntaresMediumProperties.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file MakeAntaresMediumProperties.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimMediumProperties
from icecube.clsim import I3CLSimRandomValueApplyFunction
from icecube.clsim import I3CLSimRandomValueInterpolatedDistribution
from icecube.clsim import I3CLSimRandomValueRayleighScatteringCosAngle
from icecube.clsim import I3CLSimFunctionScatLenPartic
from icecube.clsim import I3CLSimFunctionRefIndexQuanFry
from icecube.clsim import I3CLSimFunctionFromTable
from icecube.clsim import I3CLSimRandomValueMixed
from icecube.clsim import I3CLSimScalarFieldConstant
from icecube.clsim import I3CLSimVectorTransformConstant

from I3Tray import I3Units

import math
from os.path import expandvars

def GetPetzoldScatteringCosAngleDistribution():
    # petzold scattering angle distribution
    deg=math.pi/180.

    petzold_powerLawIndexBeforeFirstBin = -1.346
    petzold_data_ang = [1e-9,
                        0.100*deg,   0.126*deg,   0.158*deg,   0.200*deg,   0.251*deg,
                        0.316*deg,   0.398*deg,   0.501*deg,   0.631*deg,   0.794*deg,
                        1.000*deg,   1.259*deg,   1.585*deg,   1.995*deg,   2.512*deg,
                        3.162*deg,   3.981*deg,   5.012*deg,   6.310*deg,   7.943*deg,
                        10.000*deg,  15.000*deg,  20.000*deg,  25.000*deg,  30.000*deg,
                        35.000*deg,  40.000*deg,  45.000*deg,  50.000*deg,  55.000*deg,
                        60.000*deg,  65.000*deg,  70.000*deg,  75.000*deg,  80.000*deg,
                        85.000*deg,  90.000*deg,  95.000*deg, 100.000*deg, 105.000*deg,
                        110.000*deg, 115.000*deg, 120.000*deg, 125.000*deg, 130.000*deg,
                        135.000*deg, 140.000*deg, 145.000*deg, 150.000*deg, 155.000*deg,
                        160.000*deg, 165.000*deg, 170.000*deg, 175.000*deg, 180.000*deg]
    petzold_data_val=[0.,  # <- this 0 will be replaced right after the definition
                      1.767e+03, 1.296e+03, 9.502e+02, 6.991e+02, 5.140e+02,
                      3.764e+02, 2.763e+02, 2.188e+02, 1.444e+02, 1.022e+02,
                      7.161e+01, 4.958e+01, 3.395e+01, 2.281e+01, 1.516e+01,
                      1.002e+01, 6.580e+00, 4.295e+00, 2.807e+00, 1.819e+00,
                      1.153e+00, 4.893e-01, 2.444e-01, 1.472e-01, 8.609e-02,
                      5.931e-02, 4.210e-02, 3.067e-02, 2.275e-02, 1.699e-02,
                      1.313e-02, 1.046e-02, 8.488e-03, 6.976e-03, 5.842e-03,
                      4.953e-03, 4.292e-03, 3.782e-03, 3.404e-03, 3.116e-03,
                      2.912e-03, 2.797e-03, 2.686e-03, 2.571e-03, 2.476e-03,
                      2.377e-03, 2.329e-03, 2.313e-03, 2.365e-03, 2.506e-03,
                      2.662e-03, 2.835e-03, 3.031e-03, 3.092e-03, 3.154e-03]
    # power law for values below bin 1
    petzold_data_val[0] = 2.*math.pi*math.sin(petzold_data_ang[1])*petzold_data_val[1] * ((petzold_data_ang[0]/petzold_data_ang[1])**petzold_powerLawIndexBeforeFirstBin)
    
    for i in range(len(petzold_data_val)):
        petzold_data_val[i] = 2.*math.pi*petzold_data_val[i]*math.sin(petzold_data_ang[i])
    
    #petzoldScatModelCosAngle = I3CLSimRandomValueTabulatedDistributionCosAngle(angles=petzold_data_ang[1:], values=petzold_data_val[1:], powerLawIndexBeforeFirstBin=petzold_powerLawIndexBeforeFirstBin)

    petzoldScatModelAngle = I3CLSimRandomValueInterpolatedDistribution(x=petzold_data_ang, y=petzold_data_val)
    petzoldScatModelCosAngle = I3CLSimRandomValueApplyFunction(petzoldScatModelAngle, "cos")
    
    return petzoldScatModelCosAngle

def GetAntaresScatteringCosAngleDistribution():
    petzoldScatModel = GetPetzoldScatteringCosAngleDistribution()
    rayleighScatModel = I3CLSimRandomValueRayleighScatteringCosAngle()
    antaresScatModel = I3CLSimRandomValueMixed(
        firstDistribution=rayleighScatModel,
        secondDistribution=petzoldScatModel,
        fractionOfFirstDistribution=0.17)
    return antaresScatModel


def MakeAntaresMediumProperties():
    m = I3CLSimMediumProperties(mediumDensity=1.039*I3Units.g/I3Units.cm3,
                                layersNum=1, 
                                layersZStart=-310.*I3Units.m, 
                                layersHeight=2500.*I3Units.m,
                                rockZCoordinate=-310.*I3Units.m,
                                airZCoordinate=2500.*I3Units.m-310.*I3Units.m)
    
    antaresScatModel = GetAntaresScatteringCosAngleDistribution()
    m.SetScatteringCosAngleDistribution(antaresScatModel)

    # no ice/water anisotropy. all of these three are no-ops
    m.SetDirectionalAbsorptionLengthCorrection(I3CLSimScalarFieldConstant(1.))
    m.SetPreScatterDirectionTransform(I3CLSimVectorTransformConstant())
    m.SetPostScatterDirectionTransform(I3CLSimVectorTransformConstant())

    # no "ice" tilt
    m.SetIceTiltZShift(I3CLSimScalarFieldConstant(0.))

    antaresScattering = I3CLSimFunctionScatLenPartic(volumeConcentrationSmallParticles=0.0075*I3Units.perMillion, volumeConcentrationLargeParticles=0.0075*I3Units.perMillion)
    antaresPhaseRefIndex = I3CLSimFunctionRefIndexQuanFry(pressure=215.82225*I3Units.bar, temperature=13.1, salinity=38.44*I3Units.perThousand)
    
    # these are for ANTARES (mix of Smith&Baker water and Antares site measurements)
    # absorption lengths starting from 290nm in 10nm increments
    absLenTable = [  4.65116279,   7.1942446,    9.17431193,  10.57082452,  12.62626263,  14.08450704,
                    15.89825119,  18.93939394,  21.14164905,  24.09638554,  27.54820937,
                    30.76923077,  34.36426117,  39.21568627,  42.19409283,  45.87155963,  50.,
                    52.35602094,  54.94505495,  54.94505495,  51.02040816,  38.91050584,
                    28.01120448,  20.96436059,  19.72386588,  17.92114695,  15.67398119,
                    14.12429379,  12.51564456,   9.25925926,   6.36942675,   4.09836066,
                     3.46020761]
    antaresAbsorption = I3CLSimFunctionFromTable(290.*I3Units.nanometer, 10.*I3Units.nanometer, absLenTable)
    
    m.SetAbsorptionLength(0, antaresAbsorption)
    m.SetScatteringLength(0, antaresScattering)
    m.SetPhaseRefractiveIndex(0, antaresPhaseRefIndex)
    
    return m

