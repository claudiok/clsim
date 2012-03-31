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

# do not build clsim if there is no OpenCL on this system
if(OPENCL_FOUND)

  # the IceTray project definition
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
    private/clsim/I3CLSimLightSourceParameterization.cxx
    private/clsim/I3CLSimLightSourceToStepConverterPPC.cxx
    private/clsim/I3CLSimLightSourceToStepConverter.cxx
    private/clsim/I3CLSimLightSourceToStepConverterUtils.cxx
    private/clsim/I3CLSimPhoton.cxx
    private/clsim/I3CLSimPhotonHistory.cxx
    private/clsim/I3CLSimPMTPhotonSimulator.cxx
    private/clsim/random_value/I3CLSimRandomValueApplyFunction.cxx
    private/clsim/random_value/I3CLSimRandomValue.cxx
    private/clsim/random_value/I3CLSimRandomValueHenyeyGreenstein.cxx
    private/clsim/random_value/I3CLSimRandomValueInterpolatedDistribution.cxx
    private/clsim/random_value/I3CLSimRandomValueMixed.cxx
    private/clsim/random_value/I3CLSimRandomValueRayleighScatteringCosAngle.cxx
    private/clsim/random_value/I3CLSimRandomValueSimplifiedLiu.cxx
    private/clsim/random_value/I3CLSimRandomValueWlenCherenkovNoDispersion.cxx
    private/clsim/I3CLSimSimpleGeometryFromI3Geometry.cxx
    private/clsim/I3CLSimSimpleGeometryTextFile.cxx
    private/clsim/I3CLSimSimpleGeometryUserConfigurable.cxx
    private/clsim/I3CLSimStep.cxx
    private/clsim/function/I3CLSimFunctionAbsLenIceCube.cxx
    private/clsim/function/I3CLSimFunctionConstant.cxx
    private/clsim/function/I3CLSimFunction.cxx
    private/clsim/function/I3CLSimFunctionFromTable.cxx
    private/clsim/function/I3CLSimFunctionPolynomial.cxx
    private/clsim/function/I3CLSimFunctionRefIndexIceCube.cxx
    private/clsim/function/I3CLSimFunctionRefIndexQuanFry.cxx
    private/clsim/function/I3CLSimFunctionScatLenIceCube.cxx
    private/clsim/function/I3CLSimFunctionScatLenPartic.cxx
    private/clsim/I3CLSimLightSource.cxx
    private/clsim/I3CLSimSpectrumTable.cxx
    private/clsim/I3CLSimFlasherPulse.cxx
    private/clsim/I3CLSimLightSourceToStepConverterFlasher.cxx
    private/clsim/util/I3MuonSlicer.cxx
    private/clsim/I3Photon.cxx
    private/clsim/dom/I3PhotonToMCHitConverter.cxx
    private/clsim/shadow/I3ShadowedPhotonRemover.cxx
    private/clsim/shadow/I3ShadowedPhotonRemoverModule.cxx
    private/clsim/shadow/I3ExtraGeometryItem.cxx
    private/clsim/shadow/I3ExtraGeometryItemUnion.cxx
    private/clsim/shadow/I3ExtraGeometryItemMove.cxx
    private/clsim/shadow/I3ExtraGeometryItemCylinder.cxx

    # tableio converters
    private/clsim/converter/I3PhotonConverter.cxx
    private/clsim/converter/I3MCHitConverterWithIDs.cxx

    # private/geant4
    private/geant4/I3CLSimLightSourceToStepConverterGeant4.cxx
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
    private/test/I3CLSimFunctionTester.cxx
  )
  SET(LIB_${PROJECT_NAME}_PROJECTS icetray dataio dataclasses phys-services tableio simclasses)
  SET(LIB_${PROJECT_NAME}_TOOLS python boost log4cplus gsl clhep geant4 opencl)

  # The icecube master header file I3.h seems to include inttypes.h without setting
  # __STDC_FORMAT_MACROS. clsim makes use of the PRI* preprocessor definitions,
  # so it needs those macros. Add the definition right here to work around this.
  ADD_DEFINITIONS(-D__STDC_FORMAT_MACROS)

  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/clsim/resources/safeprimes_base32.txt)
    colormsg(RED   "+-- Downloading and unpacking the safeprimes_base32.txt data file...")
    colormsg(RED   "    (this file will only be downloaded once)")
    EXECUTE_PROCESS(COMMAND ${CMAKE_SOURCE_DIR}/clsim/resources/download_safeprimes.sh ${CMAKE_SOURCE_DIR}/clsim/resources)
    colormsg(GREEN "+-- ...download complete.")
  else(NOT EXISTS ${CMAKE_SOURCE_DIR}/clsim/resources/safeprimes_base32.txt)
    colormsg(CYAN  "+-- safeprimes_base32.txt already downloaded")
  endif(NOT EXISTS ${CMAKE_SOURCE_DIR}/clsim/resources/safeprimes_base32.txt)

  # check for Q-Frame support
  if(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)
    # we assume that the existence of this module means that we have
    # to deal with a Q-frame-enabled meta-project
    ADD_DEFINITIONS(-DIS_Q_FRAME_ENABLED)
    colormsg(GREEN "+-- Q-frame-enabled")
  else(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)
    colormsg(YELLOW "+-- old-style, non-Q-frame module")
  endif(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)


  # add some things if hit-maker is available (IceCube-specific)
  if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/hit-maker)
    LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS hit-maker)
    #LIST(APPEND LIB_${PROJECT_NAME}_TOOLS nothing) # no extra tools
    LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
         private/clsim/I3CLSimPMTPhotonSimulatorIceCube.cxx
         )

    ADD_DEFINITIONS(-DUSE_HIT_MAKER)
    colormsg(GREEN "+-- hit-maker support (for pre-Rosencrantz simulations)")
  endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/hit-maker)


  # heck for the existence of icetray/python/dataclass_suite.hpp
  # If this file does not exist, a local copy of the python
  # binding helper headers is used. TODO: This needs to go away
  # as soon as an offline-software release contains these headers.
  # (They already are in icetray/trunk.)
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/python/dataclass_suite.hpp)
    ADD_DEFINITIONS(-DNO_PYTHON_DATACLASS_SUITE)
    colormsg(YELLOW "+-- icetray does not have dataclass_suite.hpp. Using local version.")
  endif(NOT EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/python/dataclass_suite.hpp)


  # check if multiPMT support is available and add multiPMT-only 
  # source files (this is still for the KM3NeT multiPMT version of IceTray)
  # support for 3-item OMKeys is coming soon!
  if(EXISTS ${CMAKE_SOURCE_DIR}/dataclasses/public/dataclasses/geometry/I3OMTypeInfo.h)
    LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
         private/clsim/dom/I3PhotonToMCHitConverterForMultiPMT.cxx
         )

    ADD_DEFINITIONS(-DHAS_MULTIPMT_SUPPORT)
    colormsg(GREEN "+-- KM3NeT old-style multiPMT support")
  endif(EXISTS ${CMAKE_SOURCE_DIR}/dataclasses/public/dataclasses/geometry/I3OMTypeInfo.h)

  # check for glshovel support
  if(QT4_FOUND AND NOT QT_CONFIG_ERROR AND OPENGL_FOUND AND GLUT_FOUND)
    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/glshovel)
      # We seem to have glshovel and its requirements available.
      # Build the I3Photon renderer!
      LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS glshovel)
      LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
           private/glshovel/render/PhotonPaths.cxx
       )
      include(${QT_USE_FILE})

      colormsg(GREEN "+-- glshovel I3Photon renderers")
    endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/glshovel)
  endif(QT4_FOUND AND NOT QT_CONFIG_ERROR AND OPENGL_FOUND AND GLUT_FOUND)


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
    colormsg(GREEN "+-- gmp support (make_safeprimes utility)")
  else(GMP_FOUND)
    colormsg(CYAN  "+-- no gmp support (make_safeprimes utility)")
  endif(GMP_FOUND)

  add_subdirectory(private/pybindings)

else(OPENCL_FOUND) # there is no OpenCL installed on this system
    colormsg(RED   "+-- OpenCL is not installed on your system. Cannot build clsim!")
endif(OPENCL_FOUND)
