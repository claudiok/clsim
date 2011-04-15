from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimMediumProperties, \
                          I3CLSimRandomValueMixed, \
                          I3CLSimRandomValueHenyeyGreenstein, \
                          I3CLSimRandomValueSimplifiedLiu, \
                          I3CLSimWlenDependentValueRefIndexIceCube, \
                          I3CLSimWlenDependentValueAbsLenIceCube, \
                          I3CLSimWlenDependentValueScatLenIceCube

from I3Tray import I3Units

import numpy, math
from os.path import expandvars


def MakeIceCubeMediumProperties(detectorCenterDepth = 1948.07*I3Units.m):
    ### read ice information from PPC-compatible tables
    
    icemodel_dat = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/ice/spice_mie/icemodel.dat"), unpack=True)
    icemodel_par = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/ice/spice_mie/icemodel.par"))
    icemodel_cfg = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/ice/spice_mie/cfg.txt"))
    
    alpha = icemodel_par[0][0]
    kappa = icemodel_par[1][0]
    A     = icemodel_par[2][0]
    B     = icemodel_par[3][0]
    D     = icemodel_par[4][0]
    E     = icemodel_par[5][0]

        
    oversizeScaling       = icemodel_cfg[0] # currently ignored
    efficiencyCorrection  = icemodel_cfg[1] # currently ignored
    liuScatteringFraction = icemodel_cfg[2]
    meanCosineTheta       = icemodel_cfg[3]

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
        print "The ice table's delta_tau values do not seem to correpsond to the equation. Loading table anyway."
    
    
    # some sanity checks on the layer files
    if len(depth)<2: raise RuntimeError("There is only a single layer in your layer definition file")
    
    layerHeight = depth[1]-depth[0]
    if layerHeight<=0.: raise RuntimeError("ice layer depths are not in increasing order")
    
    for i in range(len(depth)-1):
        if abs((depth[i+1]-depth[i]) - layerHeight) > 1e-5:
            raise RuntimeError("ice layers are not spaced evenly")
    
    
    layerZStart = []
    layerZEnd = []
    
    for z in detectorCenterDepth - depth[::-1]: # z coordinates in reverse order
        layerZStart.append(z-layerHeight/2.)
        layerZEnd.append(z+layerHeight/2.)
    
    layerZStart = numpy.array(layerZStart)
    layerZEnd = numpy.array(layerZEnd)

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
    
    phaseRefIndex = I3CLSimWlenDependentValueRefIndexIceCube(mode="phase")
    groupRefIndex = I3CLSimWlenDependentValueRefIndexIceCube(mode="group")
    for i in range(len(layerZStart)):
        m.SetPhaseRefractiveIndex(i, phaseRefIndex)

        # the IceCube group refractive index parameterization is not exactly 
        # what you would expect from the phase refractive index.
        # To use calculated values instead of the parameterization,
        # just comment this line:
        m.SetGroupRefractiveIndexOverride(i, groupRefIndex)
        
        absLen = I3CLSimWlenDependentValueAbsLenIceCube(kappa=kappa, A=A, B=B, D=D, E=E,
                                                        aDust400=a_dust400[::-1][i],   # reverse order
                                                        deltaTau=delta_tau[::-1][i])   # reverse order
        m.SetAbsorptionLength(i, absLen)

        scatLen = I3CLSimWlenDependentValueScatLenIceCube(alpha=alpha,
                                                          be400=b_e400[::-1][i])       # reverse order
        m.SetScatteringLength(i, scatLen)

    return m
