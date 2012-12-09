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
# @file CMakeLists.txt
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

if(BUILD_PYBINDINGS)
  # these files will be compiled in any case (with or without OpenCL/Geant4)
  SET(LIB_${PROJECT_NAME}_PYBINDINGS_SOURCEFILES
    I3Photon.cxx
    I3CLSimEventStatistics.cxx
    I3CLSimPMTPhotonSimulator.cxx
    I3CLSimFlasherPulse.cxx
    I3Converters.cxx
    I3ShadowedPhotonRemover.cxx
    I3ExtraGeometryItem.cxx
    module.cxx
  )

  if(BUILD_CLSIM_DATACLASSES_ONLY)
    ADD_DEFINITIONS(-DBUILD_CLSIM_DATACLASSES_ONLY)
  else(BUILD_CLSIM_DATACLASSES_ONLY)
    # add all extra source files that do depend on OpenCL and/or Geant4
    LIST(APPEND LIB_${PROJECT_NAME}_PYBINDINGS_SOURCEFILES
      I3CLSimStep.cxx
      I3CLSimPhoton.cxx
      I3CLSimPhotonHistory.cxx
      I3CLSimFunction.cxx
      I3CLSimScalarField.cxx
      I3CLSimRandomValue.cxx
      I3CLSimMediumProperties.cxx
      I3CLSimLightSourceToStepConverter.cxx
      I3CLSimStepToPhotonConverter.cxx
      I3CLSimSimpleGeometry.cxx
      I3CLSimLightSourceParameterization.cxx
      I3CLSimTester.cxx
      I3CLSimModuleHelper.cxx
      I3CLSimLightSourceToStepConverterUtils.cxx
      I3CLSimOpenCLDevice.cxx
      I3CLSimLightSource.cxx
      I3CLSimSpectrumTable.cxx
    )
  endif(BUILD_CLSIM_DATACLASSES_ONLY)

  if(BUILD_TABULATOR)
    LIST(APPEND LIB_${PROJECT_NAME}_PYBINDINGS_SOURCEFILES
      I3CLSimTabulator.cxx
    )
  endif(BUILD_TABULATOR)

  i3_add_pybindings(clsim
    ${LIB_${PROJECT_NAME}_PYBINDINGS_SOURCEFILES}
    USE_TOOLS boost python log4cplus
    USE_PROJECTS icetray dataclasses clsim
    )
endif(BUILD_PYBINDINGS)

