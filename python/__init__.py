from icecube.load_pybindings import load_pybindings
import icecube.icetray, icecube.dataclasses # be nice and pull in our dependencies
load_pybindings(__name__,__path__)

from icecube import icetray, dataclasses

from I3Tray import I3Units

import numpy, math
from os.path import expandvars


# petzold scattering angle distribution
deg=math.pi/180.
petzold_data_ang = [0.100*deg,   0.126*deg,   0.158*deg,   0.200*deg,   0.251*deg,
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
petzold_data_val=[1.767e+03, 1.296e+03, 9.502e+02, 6.991e+02, 5.140e+02,
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
petzold_powerLawIndexBeforeFirstBin = -1.346

def MakeAntaresMediumProperties():
    m = I3CLSimMediumProperties(mediumDensity=1.039*I3Units.g/I3Units.cm3,
                                layersNum=1, 
                                layersZStart=-310.*I3Units.m, 
                                layersHeight=2500.*I3Units.m,
                                rockZCoordinate=-310.*I3Units.m,
                                airZCoordinate=2500.*I3Units.m-310.*I3Units.m)
    
    petzoldScatModel = I3CLSimRandomValueTabulatedDistributionCosAngle(angles=petzold_data_ang, values=petzold_data_val, powerLawIndexBeforeFirstBin=petzold_powerLawIndexBeforeFirstBin)
    rayleighScatModel = I3CLSimRandomValueRayleighScatteringCosAngle()
    antaresScatModel = I3CLSimRandomValueMixed(
        firstDistribution=rayleighScatModel,
        secondDistribution=petzoldScatModel,
        fractionOfFirstDistribution=0.17)
    m.SetScatteringCosAngleDistribution(antaresScatModel)
    
    antaresScattering = I3CLSimWlenDependentValueScatLenPartic(volumeConcentrationSmallParticles=0.0075*I3Units.perMillion, volumeConcentrationLargeParticles=0.0075*I3Units.perMillion)
    antaresPhaseRefIndex = I3CLSimWlenDependentValueRefIndexQuanFry(pressure=215.82225*I3Units.bar, temperature=13.1, salinity=38.44*I3Units.perThousand)
    
    # these are for ANTARES (mix of Smith&Baker water and Antares site measurements)
    # absorption lengths starting from 290nm in 10nm increments
    absLenTable = [  4.65116279,   7.1942446,    9.17431193,  10.57082452,  12.62626263,  14.08450704,
                    15.89825119,  18.93939394,  21.14164905,  24.09638554,  27.54820937,
                    30.76923077,  34.36426117,  39.21568627,  42.19409283,  45.87155963,  50.,
                    52.35602094,  54.94505495,  54.94505495,  51.02040816,  38.91050584,
                    28.01120448,  20.96436059,  19.72386588,  17.92114695,  15.67398119,
                    14.12429379,  12.51564456,   9.25925926,   6.36942675,   4.09836066,
                     3.46020761]
    antaresAbsorption = I3CLSimWlenDependentValueFromTable(290.*I3Units.nanometer, 10.*I3Units.nanometer, absLenTable)
    
    m.SetAbsorptionLength(0, antaresAbsorption)
    m.SetScatteringLength(0, antaresScattering)
    m.SetPhaseRefractiveIndex(0, antaresPhaseRefIndex)
    
    return m


def MakeIceCubeMediumProperties(detectorCenterDepth = 1948.07*I3Units.m):
    ### read ice information from PPC-compatible tables
    
    icemodel_dat = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/ice/spice_mie/icemodel.dat"), unpack=True)
    icemodel_par = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/ice/spice_mie/icemodel.par"))
    
    alpha = icemodel_par[0][0]
    kappa = icemodel_par[1][0]
    A     = icemodel_par[2][0]
    B     = icemodel_par[3][0]
    D     = icemodel_par[4][0]
    E     = icemodel_par[5][0]
    
    
    depth = icemodel_dat[0]*I3Units.m
    b_e400 = icemodel_dat[1]
    a_dust400 = icemodel_dat[2]
    delta_tau = icemodel_dat[3]
    
    # check delta_tau values against formula
    # According to the IceCube paper, this should be calculated like this.
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
        firstDistribution=I3CLSimRandomValueHenyeyGreenstein(meanCosine=0.9),
        secondDistribution=I3CLSimRandomValueSimplifiedLiu(meanCosine=0.9),
        fractionOfFirstDistribution=0.45)
    m.SetScatteringCosAngleDistribution(iceCubeScatModel)
    
    phaseRefIndex = I3CLSimWlenDependentValueRefIndexIceCube()
    for i in range(len(layerZStart)):
        m.SetPhaseRefractiveIndex(i, phaseRefIndex)

        absLen = I3CLSimWlenDependentValueAbsLenIceCube(kappa=kappa, A=A, B=B, D=D, E=E,
                                                        aDust400=a_dust400[i], 
                                                        deltaTau=delta_tau[i])
        m.SetAbsorptionLength(i, absLen)

        scatLen = I3CLSimWlenDependentValueScatLenIceCube(alpha=alpha,
                                                          be400=b_e400[i])
        m.SetScatteringLength(i, scatLen)

    return m
