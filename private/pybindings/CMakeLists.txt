#
#  Copyright (C) 2011
#  Claudio Kopper <claudio.kopper@physik.uni-erlangen.de>
#  
#  clsim is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>
#  

if(BUILD_PYBINDINGS)
  i3_add_pybindings(clsim
    I3CLSimStep.cxx
    I3CLSimPhoton.cxx
    I3CLSimPhotonHistory.cxx
    I3CLSimWlenDependentValue.cxx
    I3CLSimRandomValue.cxx
    I3CLSimMediumProperties.cxx
    I3CLSimLightSourceToStepConverter.cxx
    I3CLSimStepToPhotonConverter.cxx
    I3CLSimSimpleGeometry.cxx
    I3CLSimLightSourceParameterization.cxx
    I3CLSimTester.cxx
    I3CLSimModuleHelper.cxx
    I3CLSimLightSourceToStepConverterUtils.cxx
    I3Photon.cxx
    I3CLSimEventStatistics.cxx
    I3CLSimPMTPhotonSimulator.cxx
    I3CLSimOpenCLDevice.cxx
    I3CLSimLightSource.cxx
    I3CLSimSpectrumTable.cxx
    I3CLSimFlasherPulse.cxx
    I3Converters.cxx
    I3ShadowedPhotonRemover.cxx
    module.cxx
    USE_TOOLS boost python log4cplus
    USE_PROJECTS icetray dataclasses clsim
    )
endif(BUILD_PYBINDINGS)

