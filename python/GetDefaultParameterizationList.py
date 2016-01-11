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
# $Id: GetDefaultParameterizationList.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetDefaultParameterizationList.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimLightSourceParameterization

def GetDefaultParameterizationList(theConverter, muonOnly=False):
    fromEnergy=0.
    toEnergy=float('Inf')

    muons    = [dataclasses.I3Particle.MuMinus,
                dataclasses.I3Particle.MuPlus]

    taus     = [dataclasses.I3Particle.TauMinus,
                dataclasses.I3Particle.TauPlus]

    cascades = [dataclasses.I3Particle.Neutron,
                dataclasses.I3Particle.Hadrons,
                dataclasses.I3Particle.Pi0,
                dataclasses.I3Particle.PiPlus,
                dataclasses.I3Particle.PiMinus,
                dataclasses.I3Particle.K0_Long,
                dataclasses.I3Particle.KPlus,
                dataclasses.I3Particle.KMinus,
                dataclasses.I3Particle.PPlus,
                dataclasses.I3Particle.PMinus,
                dataclasses.I3Particle.K0_Short,
                dataclasses.I3Particle.EMinus,
                dataclasses.I3Particle.EPlus,
                dataclasses.I3Particle.Gamma,
                dataclasses.I3Particle.Brems,
                dataclasses.I3Particle.DeltaE,
                dataclasses.I3Particle.PairProd,
                dataclasses.I3Particle.NuclInt]

    parameterizationsMuon = []
    for type in muons:
        converter = \
          I3CLSimLightSourceParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=True)
        parameterizationsMuon.append(converter)

    # treat taus as muons for purposes of the converter
    for type in taus:
        converter = \
          I3CLSimLightSourceParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=True)
        parameterizationsMuon.append(converter)

    if muonOnly:
        return parameterizationsMuon

    parameterizationsOther = []
    for type in cascades:
        converter = \
          I3CLSimLightSourceParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=False)
        parameterizationsOther.append(converter)

    return parameterizationsMuon+parameterizationsOther
