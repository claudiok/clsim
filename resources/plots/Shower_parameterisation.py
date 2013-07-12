#!/usr/bin/env python

from __future__ import print_function

import math

import matplotlib
matplotlib.use("PDF")


fig_size = [8.3,11.7] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

import pylab
import scipy
import numpy

import scipy.interpolate
import scipy.integrate
import scipy.optimize

from copy import copy, deepcopy

prng = numpy.random.RandomState()


##### new shower profile (MC)

def scatterDirectionByAngle2(cost, sint, p_uhat, randomNumber):
    temp    = math.sqrt(1.0 - p_uhat[2]*p_uhat[2])

    # cosine and sine of the azimuthal angle psi.
    b = 2.0*math.pi*randomNumber
    sinp = math.sin(b)
    cosp = math.cos(b)
    
    
    if temp == 0.0:
      # normal incident.
      p_uhat_x = sint*cosp
      p_uhat_y = sint*sinp
      p_uhat_z = math.copysign(cost,p_uhat[2]*cost)
    else:
      # regular incident.
      p_uhat_x = (sint*(p_uhat[0]*p_uhat[2]*cosp - p_uhat[1]*sinp))/temp + p_uhat[0]*cost
      p_uhat_y = (sint*(p_uhat[1]*p_uhat[2]*cosp + p_uhat[0]*sinp))/temp + p_uhat[1]*cost
      p_uhat_z = -sint*cosp*temp + p_uhat[2]*cost
    
    return [p_uhat_x, p_uhat_y, p_uhat_z]

def scatterDirectionByAngle(cosa, sina, direction, randomNumber):
    # randomize direction of scattering (rotation around old direction axis)
    b=2.0*math.pi*randomNumber
    cosb=numpy.cos(b)
    sinb=numpy.sin(b)
    
    # Rotate new direction into absolute frame of reference 
    sinth = numpy.sqrt(max(0.0, 1.0-(direction[2])**2.))
    
    newDirection=[None,None,None]
    
    if sinth>0.:  # Current direction not vertical, so rotate 
        newDirection[0]=direction[0]*cosa-(sina*(direction[1]*cosb+direction[2]*direction[0]*sinb))/sinth
        newDirection[1]=direction[1]*cosa+(sina*(direction[0]*cosb-direction[2]*direction[1]*sinb))/sinth
        newDirection[2]=direction[2]*cosa+(sina*sinb*sinth);
    else:         # Current direction is vertical, so this is trivial
        newDirection[0]=sina*cosb
        newDirection[1]=sina*sinb
        newDirection[2]=cosa*numpy.sign(direction[2])

    #recip_length = 1./math.sqrt((newDirection[0]**2.) + (newDirection[1]**2.) + (newDirection[2]**2.))
    #newDirection[0] *= recip_length
    #newDirection[1] *= recip_length
    #newDirection[2] *= recip_length
    
    return newDirection

def getPhaseRefIndex(wavelength):
    # Quan&Fry (taken from W. Schuster's thesis):
    refind_S  = 38.44    # salinity in ppt
    refind_T  = 13.1     # temperature in degC
    refind_P  = 213.0    # ambient pressure [atm]   # 213 bar in comb. with the previous salinity and temp. seems to approximate the km3 tables very closely

    refind_n0 = 1.31405  # offset
    refind_n1 = 1.45e-5
    refind_n2 = 1.779e-4
    refind_n3 = 1.05e-6
    refind_n4 = 1.6e-8
    refind_n5 = 2.02e-6
    refind_n6 = 15.868
    refind_n7 = 0.01155
    refind_n8 = 0.00423
    refind_n9 = 4382.
    refind_n10 = 1.1455e6
    # these get used in the calculation:
    refind_a0 = refind_n0+(refind_n2-refind_n3*refind_T+refind_n4*refind_T*refind_T)*refind_S-refind_n5*refind_T*refind_T
    refind_a1 = refind_n1
    refind_a2 = refind_n6+refind_n7*refind_S-refind_n8*refind_T
    refind_a3 = -refind_n9
    refind_a4 = refind_n10

    x = 1./wavelength
    return (refind_a0  +  refind_a1*refind_P  +  x*(refind_a2 + x*(refind_a3 + x*refind_a4)))


#
def makeShowerParticleDirection(photonDirection, a=0.39, b=2.61):
    #const float a=0.39, b=2.61;
    #const float I=1-expf(-b*powf(2, a));
    #float cs=max(1-powf(-logf(1-xrnd(s)*I)/b, 1/a), -1.0f);
    #float si=sqrtf(1-cs*cs); rotate(cs, si, n, s);

    # rotate to shower particle direction
    I = 1.0 - numpy.exp(-b*(2.**a))
    cs = max(1.0-((-numpy.log(1.-prng.uniform(0.,1.)*I)/b)**(1./a)), -1.0)
    si = numpy.sqrt(1.0-cs**2)
    photonDirection1 = scatterDirectionByAngle(cs, si, photonDirection, prng.uniform(0.,1.))

    return photonDirection1

def makeShoweParticleAngle(photonDirection, a=0.39, b=2.61):
    photonDirection2 = makeShowerParticleDirection(photonDirection, a, b)
    return photonDirection2[0]*photonDirection[0] + photonDirection2[1]*photonDirection[1] + photonDirection2[2]*photonDirection[2]

def makeShowerPhotonDirection(photonDirection, a=0.39, b=2.61):
    #photonDirection = [0.,0.,1.] # initial direction

    MEDIUM_WLEN_NUM=17
    MEDIUM_WLEN_START=290.
    MEDIUM_WLEN_STEP=20.
    MEDIUM_MIN_PHOTON_ENERGY=1./float(MEDIUM_WLEN_START+MEDIUM_WLEN_NUM*MEDIUM_WLEN_STEP)
    MEDIUM_MAX_PHOTON_ENERGY=1./float(MEDIUM_WLEN_START)

    # make a photon energy
    photonEnergy = MEDIUM_MIN_PHOTON_ENERGY + prng.uniform(0.,1.) * (MEDIUM_MAX_PHOTON_ENERGY-MEDIUM_MIN_PHOTON_ENERGY)
    photonWavelength = 1./photonEnergy
    
    # find the phase refractive index
    cosCherenkov = 1./getPhaseRefIndex(photonWavelength)
    sinCherenkov = math.sqrt(1.0-(cosCherenkov**2.));
    
    
    # fixed cherenkov angle for now
    #cosCherenkov = 1./1.3499 # Antares "official" software value
    #cosCherenkov = 1./1.348988618 # g4sim value @ 470nm
    cosCherenkov = 0.74
    sinCherenkov = math.sqrt(1.-cosCherenkov**2.)
    
    photonDirection1 = makeShowerParticleDirection(photonDirection, a, b)

    # and now rotate to cherenkov emission direction
    photonDirection2 = scatterDirectionByAngle(cosCherenkov, sinCherenkov, photonDirection1, prng.uniform(0.,1.))

    #print photonDirection[0]

    return photonDirection2[0]*photonDirection[0] + photonDirection2[1]*photonDirection[1] + photonDirection2[2]*photonDirection[2]


def genMCHistograms(generator, samples=100000, numBins=1000):
    cosAngles = []
    for i in range(samples):
        cosAngles.append(generator())
    cosAngles = numpy.array(cosAngles) # convert to numpy array
    numCos_orig, binsCos = scipy.histogram(cosAngles, range=(-1.,1.), bins=numBins)
    del cosAngles # not needed anymore

    numCos=[]
    for i, number in enumerate(numCos_orig):
        numCos.append(float(number)/float(samples)/float(2./float(numBins)))
    numCos=numpy.array(numCos)

    binsCos = numpy.array(binsCos[:-1])+(binsCos[1]-binsCos[0])/2.

    return dict(num=numCos, bins=binsCos)

def showerParticleDirectionDistribution(cosAngle, a=0.39, b=2.61):
    if cosAngle < -1.: return float('NaN')
    if cosAngle > 1.: return float('NaN')
    
    a = 0.39
    b = 2.61
    I = 1.0 - numpy.exp(-b*(2.**a))
    return (1.-numpy.exp(-b*((1.-cosAngle)**a)))/I

def showerParticleDirectionDistributionDerivative(cosAngle, a=0.39, b=2.61):
    epsilon = 1e-5
    return (showerParticleDirectionDistribution(cosAngle, a, b)-showerParticleDirectionDistribution(cosAngle+epsilon, a, b))/epsilon
showerParticleDirectionDistributionDerivative = numpy.vectorize(showerParticleDirectionDistributionDerivative)

def showerParticleDirectionDistributionDerivativeAnalytical(cosAngle, a=0.39, b=2.61):
    #a = 0.39
    #b = 2.61
    #a = 0.45
    #b = 3.3
    I = 1.0 - numpy.exp(-b*(2.**a))
    return (b*numpy.exp(-b*((1.-cosAngle)**a)) * (a*(1.-cosAngle)**(a-1.)) )/I





##### older shower profiles

def xangelec(cosAngle):
    """
    Ported from geasim/shopar.f
    
    density of photon fraction photon (per sr)
    emitted in direction theta / w.r.t electron direction
    ct (cos(theta))
    The data are from 100 GeV electron showers
    for cos(th) < 0.4 a parametrisation has been used
    for cos(th) > 0.4 the original histogram is used
    with bin width 0.001 in cos(th)
    """
    
    p1=-3.75512
    p2=2.0295
    p3=0.8671
    p4=0.4939

    xcos = [
    .63240E-01,.63499E-01,.63689E-01,.64324E-01,.64445E-01,
    .64766E-01,.65007E-01,.65276E-01,.65482E-01,.65551E-01,
    .65702E-01,.65692E-01,.65632E-01,.65759E-01,.66151E-01,
    .66326E-01,.66959E-01,.67267E-01,.67309E-01,.67326E-01,
    .67611E-01,.67787E-01,.67606E-01,.68210E-01,.68313E-01,
    .68874E-01,.69105E-01,.69680E-01,.69797E-01,.70012E-01,
    .70007E-01,.70023E-01,.70225E-01,.70433E-01,.70868E-01,
    .71245E-01,.71151E-01,.71548E-01,.72021E-01,.71784E-01,
    .72057E-01,.72311E-01,.73144E-01,.73479E-01,.73933E-01,
    .74202E-01,.74397E-01,.74372E-01,.74438E-01,.74663E-01,
    .75164E-01,.75566E-01,.75925E-01,.76160E-01,.76490E-01,
    .76862E-01,.76864E-01,.77213E-01,.77450E-01,.77774E-01,
    .77944E-01,.78127E-01,.78597E-01,.79231E-01,.79459E-01,
    .79450E-01,.79501E-01,.79833E-01,.80291E-01,.80330E-01,
    .81147E-01,.81717E-01,.81856E-01,.81989E-01,.82586E-01,
    .82821E-01,.82883E-01,.83794E-01,.83963E-01,.84446E-01,
    .84847E-01,.84876E-01,.84914E-01,.84844E-01,.85334E-01,
    .85753E-01,.85869E-01,.86403E-01,.87502E-01,.87547E-01,
    .87895E-01,.88062E-01,.88534E-01,.89359E-01,.89366E-01,
    .89566E-01,.89872E-01,.90530E-01,.91211E-01,.91384E-01,
    .91996E-01,.92409E-01,.92649E-01,.92863E-01,.92748E-01,
    .93650E-01,.93550E-01,.94521E-01,.94892E-01,.95806E-01,
    .95950E-01,.96401E-01,.96830E-01,.97229E-01,.97803E-01,
    .97897E-01,.98161E-01,.98819E-01,.99659E-01,.10009E+00,
    .10016E+00,.99950E-01,.10075E+00,.10199E+00,.10195E+00,
    .10208E+00,.10263E+00,.10306E+00,.10330E+00,.10386E+00,
    .10468E+00,.10593E+00,.10625E+00,.10646E+00,.10719E+00,
    .10776E+00,.10803E+00,.10890E+00,.10907E+00,.10849E+00,
    .10936E+00,.10975E+00,.11074E+00,.11172E+00,.11226E+00,
    .11245E+00,.11278E+00,.11310E+00,.11322E+00,.11377E+00,
    .11417E+00,.11528E+00,.11643E+00,.11613E+00,.11711E+00,
    .11756E+00,.11818E+00,.11957E+00,.11974E+00,.11992E+00,
    .12013E+00,.12161E+00,.12245E+00,.12330E+00,.12363E+00,
    .12387E+00,.12389E+00,.12443E+00,.12534E+00,.12635E+00,
    .12695E+00,.12776E+00,.12810E+00,.12917E+00,.12989E+00,
    .13042E+00,.13103E+00,.13129E+00,.13258E+00,.13239E+00,
    .13356E+00,.13447E+00,.13520E+00,.13596E+00,.13656E+00,
    .13784E+00,.13915E+00,.14004E+00,.14073E+00,.14152E+00,
    .14214E+00,.14276E+00,.14348E+00,.14387E+00,.14556E+00,
    .14669E+00,.14753E+00,.14879E+00,.14887E+00,.14979E+00,
    .15065E+00,.15188E+00,.15244E+00,.15374E+00,.15454E+00,
    .15542E+00,.15627E+00,.15775E+00,.15843E+00,.16032E+00,
    .16060E+00,.16156E+00,.16360E+00,.16423E+00,.16591E+00,
    .16517E+00,.16776E+00,.16830E+00,.16867E+00,.17021E+00,
    .17053E+00,.17251E+00,.17312E+00,.17397E+00,.17603E+00,
    .17807E+00,.17935E+00,.17988E+00,.18204E+00,.18318E+00,
    .18417E+00,.18466E+00,.18617E+00,.18794E+00,.18915E+00,
    .19031E+00,.19254E+00,.19382E+00,.19654E+00,.19776E+00,
    .19881E+00,.19995E+00,.20144E+00,.20251E+00,.20449E+00,
    .20595E+00,.20873E+00,.20964E+00,.21134E+00,.21403E+00,
    .21645E+00,.21847E+00,.21995E+00,.22201E+00,.22294E+00,
    .22539E+00,.22670E+00,.22869E+00,.23105E+00,.23392E+00,
    .23557E+00,.23745E+00,.23904E+00,.24182E+00,.24351E+00,
    .24459E+00,.24757E+00,.24906E+00,.25237E+00,.25540E+00,
    .25843E+00,.26117E+00,.26243E+00,.26605E+00,.26654E+00,
    .27129E+00,.27407E+00,.27818E+00,.28156E+00,.28397E+00,
    .28857E+00,.29159E+00,.29340E+00,.29748E+00,.29953E+00,
    .30294E+00,.30702E+00,.31058E+00,.31424E+00,.31842E+00,
    .32293E+00,.32774E+00,.33076E+00,.33584E+00,.34251E+00,
    .34708E+00,.35148E+00,.35760E+00,.36300E+00,.36657E+00,
    .37252E+00,.37822E+00,.38385E+00,.39070E+00,.39555E+00,
    .40247E+00,.40876E+00,.41415E+00,.42073E+00,.42916E+00,
    .43834E+00,.44819E+00,.45648E+00,.46789E+00,.47832E+00,
    .48813E+00,.50003E+00,.51267E+00,.52784E+00,.53624E+00,
    .55019E+00,.56289E+00,.57777E+00,.58864E+00,.60589E+00,
    .62372E+00,.64105E+00,.66195E+00,.68502E+00,.70688E+00,
    .73705E+00,.77001E+00,.80268E+00,.84599E+00,.88811E+00,
    .95200E+00,.10159E+01,.11088E+01,.12279E+01,.14213E+01,
    .18671E+01,.16344E+01,.13032E+01,.11645E+01,.10600E+01,
    .98510E+00,.92322E+00,.87016E+00,.82707E+00,.79339E+00,
    .75763E+00,.72645E+00,.70639E+00,.67949E+00,.65760E+00,
    .64205E+00,.62210E+00,.60486E+00,.59491E+00,.57834E+00,
    .56501E+00,.55239E+00,.54288E+00,.52553E+00,.51559E+00,
    .50349E+00,.49247E+00,.48069E+00,.47166E+00,.46091E+00,
    .45086E+00,.44160E+00,.43619E+00,.43020E+00,.42273E+00,
    .41621E+00,.41017E+00,.40419E+00,.39610E+00,.39271E+00,
    .38698E+00,.37943E+00,.37582E+00,.36911E+00,.36422E+00,
    .35712E+00,.35421E+00,.34890E+00,.34306E+00,.33874E+00,
    .33530E+00,.33074E+00,.32811E+00,.32320E+00,.32247E+00,
    .31785E+00,.31420E+00,.30983E+00,.30670E+00,.30136E+00,
    .29736E+00,.29581E+00,.29293E+00,.28994E+00,.28710E+00,
    .28254E+00,.28070E+00,.27675E+00,.27510E+00,.27374E+00,
    .27019E+00,.26803E+00,.26647E+00,.26247E+00,.25974E+00,
    .25931E+00,.25569E+00,.25351E+00,.25230E+00,.25041E+00,
    .24593E+00,.24408E+00,.24244E+00,.23958E+00,.23677E+00,
    .23602E+00,.23394E+00,.23239E+00,.23141E+00,.22858E+00,
    .22528E+00,.22336E+00,.22208E+00,.22001E+00,.21921E+00,
    .21730E+00,.21532E+00,.21534E+00,.21128E+00,.20909E+00,
    .20873E+00,.20718E+00,.20474E+00,.20433E+00,.20256E+00,
    .20260E+00,.20041E+00,.19827E+00,.19753E+00,.19594E+00,
    .19437E+00,.19254E+00,.19160E+00,.19024E+00,.18917E+00,
    .18758E+00,.18706E+00,.18508E+00,.18372E+00,.18196E+00,
    .18100E+00,.17999E+00,.17986E+00,.17784E+00,.17590E+00,
    .17483E+00,.17394E+00,.17252E+00,.17129E+00,.17035E+00,
    .16977E+00,.16888E+00,.16763E+00,.16677E+00,.16516E+00,
    .16420E+00,.16406E+00,.16304E+00,.16157E+00,.16071E+00,
    .16046E+00,.15950E+00,.15737E+00,.15743E+00,.15654E+00,
    .15444E+00,.15449E+00,.15374E+00,.15358E+00,.15170E+00,
    .15033E+00,.15079E+00,.14952E+00,.14922E+00,.14814E+00,
    .14677E+00,.14565E+00,.14463E+00,.14478E+00,.14373E+00,
    .14244E+00,.14239E+00,.14122E+00,.14091E+00,.13955E+00,
    .13886E+00,.13870E+00,.13686E+00,.13591E+00,.13619E+00,
    .13508E+00,.13415E+00,.13385E+00,.13271E+00,.13186E+00,
    .13151E+00,.13005E+00,.12948E+00,.12871E+00,.12876E+00,
    .12830E+00,.12734E+00,.12625E+00,.12628E+00,.12538E+00,
    .12416E+00,.12402E+00,.12326E+00,.12286E+00,.12213E+00,
    .12170E+00,.12095E+00,.12025E+00,.11949E+00,.11939E+00,
    .11889E+00,.11710E+00,.11748E+00,.11641E+00,.11617E+00,
    .11533E+00,.11572E+00,.11441E+00,.11379E+00,.11350E+00,
    .11339E+00,.11261E+00,.11207E+00,.11271E+00,.11215E+00,
    .11066E+00,.11034E+00,.10964E+00,.10948E+00,.10934E+00,
    .10899E+00,.10811E+00,.10791E+00,.10724E+00,.10685E+00,
    .10619E+00,.10577E+00,.10544E+00,.10491E+00,.10429E+00,
    .10333E+00,.10290E+00,.10235E+00,.10221E+00,.10129E+00,
    .10092E+00,.10064E+00,.99867E-01,.10008E+00,.99247E-01,
    .98823E-01,.98497E-01,.97372E-01,.97476E-01,.96869E-01,
    .96504E-01,.96196E-01,.95484E-01,.95535E-01,.94785E-01,
    .94270E-01,.94148E-01,.93478E-01,.92812E-01,.92415E-01,
    .92641E-01,.91932E-01,.91399E-01,.91551E-01,.91054E-01,
    .90672E-01,.89915E-01,.89259E-01,.87689E-01,.85554E-01]

    if cosAngle < 0.4:
        return math.exp(p1+p2*cosAngle+p3*(cosAngle**2.)+p4*(cosAngle**3.))
        return 0.
    else:
        i = int((cosAngle-0.4)*1000.)
        return xcos[min(i,len(xcos)-1)]
GeasimShowerProfile = numpy.vectorize(lambda x: xangelec(x))

Geant4FixedRefIndex = numpy.loadtxt("ShowerData/data_100GeV_x100.cosine.single_refindex.brunner_like.geant4.dat", unpack=True)
Geant4FixedCherenkovAngle = numpy.loadtxt("ShowerData/data_100GeV_x100.cosine.single_refindex.brunner_like.dat", unpack=True)
Geant4Realistic = numpy.loadtxt("ShowerData/data_100GeV_x100.cosine.real_refindex.dat", unpack=True)



interpolatedGeant4FixedRefIndex = scipy.interpolate.interp1d(Geant4FixedRefIndex[0]+(Geant4FixedRefIndex[0][1]-Geant4FixedRefIndex[0][0])/2., Geant4FixedRefIndex[1],bounds_error=False)
interpolatedGeant4FixedCherenkovAngle = scipy.interpolate.interp1d(Geant4FixedCherenkovAngle[0]+(Geant4FixedCherenkovAngle[0][1]-Geant4FixedCherenkovAngle[0][0])/2., Geant4FixedCherenkovAngle[1],bounds_error=False)
interpolatedGeant4Realistic = scipy.interpolate.interp1d(Geant4Realistic[0]+(Geant4Realistic[0][1]-Geant4Realistic[0][0])/2., Geant4Realistic[1],bounds_error=False)

integralGeasimShowerProfile = scipy.integrate.romberg(GeasimShowerProfile, -1., 1.)
integralGeant4FixedRefIndex = scipy.integrate.romberg(interpolatedGeant4FixedRefIndex, -0.999, 0.999)
integralGeant4FixedCherenkovAngle = scipy.integrate.romberg(interpolatedGeant4FixedCherenkovAngle, -0.999, 0.999)
integralGeant4Realistic = scipy.integrate.romberg(interpolatedGeant4Realistic, -0.999, 0.999)


GeasimShowerProfile = numpy.vectorize(lambda x: xangelec(x)/integralGeasimShowerProfile)
interpolatedGeant4FixedRefIndex = scipy.interpolate.interp1d(Geant4FixedRefIndex[0]+(Geant4FixedRefIndex[0][1]-Geant4FixedRefIndex[0][0])/2., Geant4FixedRefIndex[1]/integralGeant4FixedRefIndex,bounds_error=False)
interpolatedGeant4FixedCherenkovAngle = scipy.interpolate.interp1d(Geant4FixedCherenkovAngle[0]+(Geant4FixedCherenkovAngle[0][1]-Geant4FixedCherenkovAngle[0][0])/2., Geant4FixedCherenkovAngle[1]/integralGeant4FixedCherenkovAngle,bounds_error=False)
interpolatedGeant4Realistic = scipy.interpolate.interp1d(Geant4Realistic[0]+(Geant4Realistic[0][1]-Geant4Realistic[0][0])/2., Geant4Realistic[1]/integralGeant4Realistic,bounds_error=False)


def DirectionProfileFromPhotonProfile(cosDirection):
    cosCherenkov = 1./1.348988618 # g4sim value @ 470nm
    sinCherenkov = math.sqrt(1.-cosCherenkov**2.)
    
    def DirectionProfileFromPhotonProfileForX(cosDirection, x):
        if x==0.: return cosDirection/cosCherenkov
        sin2pix = math.sin(2.*math.pi*x)
        
        discr = cosCherenkov**2. + 4.*(sinCherenkov**2.)*(sin2pix**2.) - 4.*sinCherenkov*sin2pix*cosDirection
        if discr < 0.: 
            print("discr=", discr)
            return
        
        ret1 = (cosCherenkov + math.sqrt(discr))/(2.*sinCherenkov*sin2pix)
        ret2 = (cosCherenkov - math.sqrt(discr))/(2.*sinCherenkov*sin2pix)
        
        print("x=%f: ret1=%f, ret2=%f" % (x, ret1, ret2))
        
    for x in numpy.linspace(0.,1.,100):
        DirectionProfileFromPhotonProfileForX(cosDirection, x)
        
#       
#DirectionProfileFromPhotonProfile(1.8)

try:
    angularParticleDirectionsWithNumPhotonsWeight_hist = numpy.loadtxt("ShowerData/angularParticleDirectionsWithNumPhotonsWeight_hist.dat")
except IOError:
    

    def makeHist(cosAnglesAndWeights, numBins=2000):
        numCos_orig, binsCos = scipy.histogram(cosAnglesAndWeights[0], range=(-1.,1.), weights=cosAnglesAndWeights[1], bins=numBins)
        sumWeights = sum(cosAnglesAndWeights[1])

        numCos=[]
        for i, number in enumerate(numCos_orig):
            numCos.append(float(number)/float(sumWeights)/float(2./float(numBins)))
        numCos=numpy.array(numCos)

        binsCos = numpy.array(binsCos[:-1])+(binsCos[1]-binsCos[0])/2.

        return numpy.array([binsCos, numCos])

    histFileNames = ["ShowerData/angularParticleDirectionsWithNumPhotonsWeight01.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight02.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight03.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight04.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight05.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight06.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight07.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight08.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight09.dat",
                     "ShowerData/angularParticleDirectionsWithNumPhotonsWeight10.dat"]

    angularParticleDirectionsWithNumPhotonsWeight_hist = None
    for histFileName in histFileNames:
        angularParticleDirectionsWithNumPhotonsWeight_hist_previous = angularParticleDirectionsWithNumPhotonsWeight_hist
        angularParticleDirectionsWithNumPhotonsWeight_hist = None

        print("Loading", histFileName)
        angularParticleDirectionsWithNumPhotonsWeight = numpy.loadtxt(histFileName, unpack=True)
        angularParticleDirectionsWithNumPhotonsWeight_hist = makeHist(angularParticleDirectionsWithNumPhotonsWeight)
        del angularParticleDirectionsWithNumPhotonsWeight

        # assumes that the same bins are used
        if angularParticleDirectionsWithNumPhotonsWeight_hist_previous is not None:
            angularParticleDirectionsWithNumPhotonsWeight_hist[1] += angularParticleDirectionsWithNumPhotonsWeight_hist_previous[1]
        
    angularParticleDirectionsWithNumPhotonsWeight_hist[1] /= float(len(histFileNames))

    numpy.savetxt("ShowerData/angularParticleDirectionsWithNumPhotonsWeight_hist.dat", angularParticleDirectionsWithNumPhotonsWeight_hist)

# fit a&b
errfunc = lambda p, x, y: numpy.log(showerParticleDirectionDistributionDerivativeAnalytical(x, a=p[0], b=p[1])) - numpy.log(y)
p0 = [0.39, 2.61] # initial guess

p1, success = scipy.optimize.leastsq(errfunc, p0[:], args=(angularParticleDirectionsWithNumPhotonsWeight_hist[0], angularParticleDirectionsWithNumPhotonsWeight_hist[1]))

print("p1=", p1)
print("success=", success)

showerDist_a = p1[0]
showerDist_b = p1[1]

# best fit currently yields a=0.44341247 b=3.30908681

### do some simple MC simulation to show that throwing dice works correctly..

print("starting MC... (part 1)")
hist_showerProfileNew = genMCHistograms(lambda: makeShowerPhotonDirection([0.,0.,1.], a=showerDist_a, b=showerDist_b), samples=10000000, numBins=2000)
print("starting MC... (part 2)")
hist_showerProfileParticleNew = genMCHistograms(lambda: makeShoweParticleAngle([0.,0.,1.], a=showerDist_a, b=showerDist_b), samples=10000000, numBins=2000)
print("...MC finished")



##### preapre figure

fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)


###### plot the plot

cosAngles=numpy.linspace(-1.,1.,10000)



l, = ax.semilogy(cosAngles, GeasimShowerProfile(cosAngles), linewidth=2, color='k', label=r"\emph{geasim}")
l.set_dashes([1,1])
l, = ax.plot(cosAngles, interpolatedGeant4FixedRefIndex(cosAngles), linewidth=1.5, color='r', label=r"fixed $\theta_\mathrm{Ch}$")
l, = ax.plot(cosAngles, interpolatedGeant4FixedCherenkovAngle(cosAngles), linewidth=1.5, color='b', label=r"fixed $n$")
l, = ax.plot(cosAngles, interpolatedGeant4Realistic(cosAngles), linewidth=1.5, color='g', label=r"realistic")
l, = ax.plot(hist_showerProfileNew["bins"], hist_showerProfileNew["num"], linewidth=1.5, color='k', label=r"MC")
ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
ax.set_xlim(-1.,1.)
ax.grid(True)
ax.legend(loc='upper left')
ax.set_xlabel(r"$\cos (\alpha)$")
ax.set_ylabel(r"a.u. $\propto$ number of photons")



l, = bx.semilogy(cosAngles, GeasimShowerProfile(cosAngles), linewidth=2, color='k')
l.set_dashes([1,1])
l, = bx.plot(cosAngles, interpolatedGeant4FixedRefIndex(cosAngles), linewidth=1.5, color='r', label=r"fixed $\theta_\mathrm{Ch}$")
l, = bx.plot(cosAngles, interpolatedGeant4FixedCherenkovAngle(cosAngles), linewidth=1.5, color='b', label=r"fixed $n$")
l, = bx.plot(cosAngles, interpolatedGeant4Realistic(cosAngles), linewidth=1.5, color='g', label=r"realistic")
l, = bx.plot(hist_showerProfileNew["bins"], hist_showerProfileNew["num"], linewidth=1.5, color='k', label=r"MC")
bx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(6))
bx.set_xlim(0.715,0.765)
bx.set_ylim(4e0/1.7,3e1/1.7)
bx.grid(True)
bx.legend(loc='upper left')
bx.set_xlabel(r"$\cos (\alpha)$")
bx.set_ylabel(r"a.u. $\propto$ number of photons")

#cosAngles=numpy.linspace(-0.99,0.99,10000)
#l, = cx.semilogy(cosAngles, showerParticleDirectionDistributionDerivative(cosAngles), linewidth=1.5, color='k', label=r"empiric formula")
l, = cx.loglog(1.-hist_showerProfileParticleNew["bins"], hist_showerProfileParticleNew["num"], linewidth=1.5, color='m', label=r"MC particle directions")
l, = cx.loglog(1.-angularParticleDirectionsWithNumPhotonsWeight_hist[0], angularParticleDirectionsWithNumPhotonsWeight_hist[1], linewidth=1.5, color='b', label=r"Geant4 e- @ 100TeV")
l, = cx.loglog(1.-cosAngles, showerParticleDirectionDistributionDerivativeAnalytical(cosAngles, a=showerDist_a, b=showerDist_b), linewidth=1.5, color='k', label=r"empiric formula")
cx.grid(True)
#cx.set_xlim(-1.,1.)
#cx.set_ylim(1e-2,1e2)
cx.set_xlabel(r"$1 - \cos (\alpha)$")
cx.set_ylabel(r"a.u. $\propto$ number of photons")
cx.legend(loc='upper right')


pylab.savefig("Shower_parameterisation.pdf", transparent=True)

