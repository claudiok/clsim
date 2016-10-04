#
# Copyright (c) 2011, 2012, 2014
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# Markus Vehring <markus.vehring@icecube.wisc.edu>
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
# $Id: GetHybridParameterizationList.py 123705 2014-09-23 13:50:57Z vehring $
# 
# @file GetHybridParameterizationList.py
# @version $Revision: 123705 $
# @date $Date: 2014-09-23 09:50:57 -0400 (Tue, 23 Sep 2014) $
# @author Claudio Kopper
#

import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimLightSourceParameterization

def GetHybridParameterizationList(PPCConverter, CrossoverEnergyEM=0.1, CrossoverEnergyHadron=30.):

    fromEnergy=0.
    toEnergy=float('Inf')

    muons    = [dataclasses.I3Particle.MuMinus,
                dataclasses.I3Particle.MuPlus]

    taus     = [dataclasses.I3Particle.TauMinus,
                dataclasses.I3Particle.TauPlus]

    em_cascades = [dataclasses.I3Particle.EMinus,
                dataclasses.I3Particle.EPlus,
                dataclasses.I3Particle.Gamma,
                dataclasses.I3Particle.Brems,
                dataclasses.I3Particle.DeltaE,
                dataclasses.I3Particle.PairProd]

    hadron_cascades = [dataclasses.I3Particle.Neutron,
                dataclasses.I3Particle.Hadrons,
                dataclasses.I3Particle.PiPlus,
                dataclasses.I3Particle.PiMinus,
                dataclasses.I3Particle.Pi0, #Pi0s are parametrized as EM cascades above crossover energy, but have mass and should be propagated with hadrons 
                dataclasses.I3Particle.K0_Long,
                dataclasses.I3Particle.KPlus,
                dataclasses.I3Particle.KMinus,
                dataclasses.I3Particle.PPlus,
                dataclasses.I3Particle.PMinus,
                dataclasses.I3Particle.K0_Short, #31% decay in Pi0 and 69% in Pi+Pi- ... so it's hadronic ... ish
                dataclasses.I3Particle.NuclInt]

    parameterizationsMuon = []
    for type in muons:
        converter = \
          I3CLSimLightSourceParameterization(
            converter=PPCConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=True)
        parameterizationsMuon.append(converter)

    # Do not parametrize taus. Let GEANT4 do it's thing.

    parameterizationsEM = []
    if not CrossoverEnergyEM is None: #if no crossover is set use GEANT4
        for type in em_cascades:
            #PPC from the CrossoverEnergy to Inf. GeV
            converter = \
              I3CLSimLightSourceParameterization(
                converter=PPCConverter,
                forParticleType=type,
                fromEnergy=CrossoverEnergyEM,
                toEnergy=toEnergy, 
                needsLength=False)
            parameterizationsEM.append(converter)

    #for hadrons split everything into two different energy ranges for PPC or Geant4
    parameterizationsHadron = []
    if not CrossoverEnergyHadron is None: #if no crossover is set use GEANT4
        for type in hadron_cascades:
            #PPC from the CrossoverEnergy to Inf. GeV
            converter = \
              I3CLSimLightSourceParameterization(
                converter=PPCConverter,
                forParticleType=type,
                fromEnergy=CrossoverEnergyHadron,
                toEnergy=toEnergy,
                needsLength=False)
            parameterizationsHadron.append(converter)

    return parameterizationsMuon + parameterizationsEM + parameterizationsHadron
