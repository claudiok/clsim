#
#  $Id$
#  
#  Copyright (C) 2011
#  Claudio Kopper  <claudio.kopper@nikhef.nl>
#  
#  This file is free software; you can redistribute it and/or modify
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
i3_project(clsim
  PYTHON_DIR python
  DOCS_DIR resources/docs
  )

SET(LIB_${PROJECT_NAME}_SOURCEFILES
    
    # private/clsim
    private/clsim/I3CLSimEventStatistics.cxx
    private/clsim/I3CLSimMediumProperties.cxx
    private/clsim/I3CLSimModule.cxx
    private/clsim/I3CLSimModuleHelper.cxx
    private/clsim/I3CLSimParticleParameterization.cxx
    private/clsim/I3CLSimParticleToStepConverterPPC.cxx
    private/clsim/I3CLSimParticleToStepConverter.cxx
    private/clsim/I3CLSimParticleToStepConverterUtils.cxx
    private/clsim/I3CLSimPhoton.cxx
    private/clsim/I3CLSimPMTPhotonSimulator.cxx
    private/clsim/I3CLSimRandomValueApplyFunction.cxx
    private/clsim/I3CLSimRandomValue.cxx
    private/clsim/I3CLSimRandomValueHenyeyGreenstein.cxx
    private/clsim/I3CLSimRandomValueInterpolatedDistribution.cxx
    private/clsim/I3CLSimRandomValueMixed.cxx
    private/clsim/I3CLSimRandomValueRayleighScatteringCosAngle.cxx
    private/clsim/I3CLSimRandomValueSimplifiedLiu.cxx
    private/clsim/I3CLSimRandomValueWlenCherenkovNoDispersion.cxx
    private/clsim/I3CLSimSimpleGeometryFromI3Geometry.cxx
    private/clsim/I3CLSimSimpleGeometryTextFile.cxx
    private/clsim/I3CLSimSimpleGeometryUserConfigurable.cxx
    private/clsim/I3CLSimStep.cxx
    private/clsim/I3CLSimWlenDependentValueAbsLenIceCube.cxx
    private/clsim/I3CLSimWlenDependentValueConstant.cxx
    private/clsim/I3CLSimWlenDependentValue.cxx
    private/clsim/I3CLSimWlenDependentValueFromTable.cxx
    private/clsim/I3CLSimWlenDependentValuePolynomial.cxx
    private/clsim/I3CLSimWlenDependentValueRefIndexIceCube.cxx
    private/clsim/I3CLSimWlenDependentValueRefIndexQuanFry.cxx
    private/clsim/I3CLSimWlenDependentValueScatLenIceCube.cxx
    private/clsim/I3CLSimWlenDependentValueScatLenPartic.cxx
    private/clsim/I3MuonSlicer.cxx
    private/clsim/I3Photon.cxx
    private/clsim/I3PhotonToMCHitConverter.cxx
    private/clsim/I3HoleIceSimulator.cxx
    private/clsim/I3HoleIceSimulatorModule.cxx
    
    # tableio converters
    private/clsim/converter/I3PhotonConverter.cxx
    private/clsim/converter/I3MCHitConverterWithIDs.cxx
    
    # private/geant4
    private/geant4/I3CLSimParticleToStepConverterGeant4.cxx
    private/geant4/TrkCerenkov.cxx
    private/geant4/TrkDetectorConstruction.cxx
    private/geant4/TrkEMPhysicsUHE.cxx
    private/geant4/TrkEnergyCut.cxx
    private/geant4/TrkEventAction.cxx
    private/geant4/TrkOpticalPhysics.cxx
    private/geant4/TrkPrimaryGeneratorAction.cxx
    private/geant4/TrkStackingAction.cxx
    private/geant4/TrkUISessionToQueue.cxx
    private/geant4/TrkUserEventInformation.cxx

    # private/opencl/
    private/opencl/I3CLSimHelperGenerateGeometrySource.cxx
    private/opencl/I3CLSimHelperGenerateMediumPropertiesSource.cxx
    private/opencl/I3CLSimHelperGenerateMediumPropertiesSource_Optimizers.cxx
    private/opencl/I3CLSimStepToPhotonConverterOpenCL.cxx
    private/opencl/I3CLSimOpenCLDevice.cxx

    # private/test/
    private/test/I3CLSimMediumPropertiesTester.cxx
    private/test/I3CLSimRandomDistributionTester.cxx
    private/test/I3CLSimTesterBase.cxx
    private/test/I3CLSimWlenDependentValueTester.cxx
    )
SET(LIB_${PROJECT_NAME}_PROJECTS icetray dataio dataclasses phys-services tableio simclasses)
SET(LIB_${PROJECT_NAME}_TOOLS python boost log4cplus gsl clhep geant4 opencl)

# The icecube master header file I3.h seems to include inttypes.h without setting
# __STDC_FORMAT_MACROS. clsim makes use of the PRI* preprocessor definitions,
# so it needs those macros. Add the definition right here to work around this.
ADD_DEFINITIONS(-D__STDC_FORMAT_MACROS)

# check for Q-Frame support
if(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)
    # we assume that the existence of this module means that we have
    # to deal with a Q-frame-enabled meta-project
    ADD_DEFINITIONS(-DIS_Q_FRAME_ENABLED)
    MESSAGE(STATUS "  + Q-frame-enabled")
else(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)
    MESSAGE(STATUS "  - old-style, non-Q-frame module")
endif(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)


# add some things if hit-maker is available (IceCube-specific)
if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/hit-maker)
    LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS hit-maker)
    #LIST(APPEND LIB_${PROJECT_NAME}_TOOLS nothing) # no extra tools
    LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
         private/clsim/I3CLSimPMTPhotonSimulatorIceCube.cxx
         )
    
    ADD_DEFINITIONS(-DUSE_HIT_MAKER)
    MESSAGE(STATUS "  + hit-maker support")
endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/hit-maker)


# check if multiPMT support is available and add multiPMT-only 
# source files
if(EXISTS ${CMAKE_SOURCE_DIR}/dataclasses/public/dataclasses/geometry/I3OMTypeInfo.h)
  LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
       private/clsim/I3PhotonToMCHitConverterForMultiPMT.cxx
       )

    ADD_DEFINITIONS(-DHAS_MULTIPMT_SUPPORT)
    MESSAGE(STATUS "  + multiPMT support")
endif(EXISTS ${CMAKE_SOURCE_DIR}/dataclasses/public/dataclasses/geometry/I3OMTypeInfo.h)


i3_add_library(${PROJECT_NAME}
    ${LIB_${PROJECT_NAME}_SOURCEFILES}
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
    )


# the make-safeprimes tool needs gmp, so only compile it if that tool is available

# add some things if hit-maker is available (IceCube-specific)
if(GMP_FOUND)
    i3_executable(make_safeprimes
      private/make_safeprimes/main.c
      USE_TOOLS python boost gmp
      )
    MESSAGE(STATUS "  + gmp support (make_safeprimes utility)")
else(GMP_FOUND)
    MESSAGE(STATUS "  - no gmp support (make_safeprimes utility)")
endif(GMP_FOUND)

add_subdirectory(private/pybindings)

