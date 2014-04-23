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

# make this file compatible with old versions of icetray that do not have "colormsg"
if(NOT COMMAND colormsg)
  function (colormsg)
    if($ARGC EQUAL 1)
      message(STATUS ${ARGV0})
    else($ARGC EQUAL 1)
      message(STATUS "  " ${ARGV1})
    endif($ARGC EQUAL 1)
  endfunction()
endif(NOT COMMAND colormsg)

# Do we have OpenCL and Geant4 installed on this system? If not,
# build a simple version of clsim with dataclasses only.
if (OPENCL_FOUND)
  set (BUILD_CLSIM_DATACLASSES_ONLY NO)
else (OPENCL_FOUND)
  set (BUILD_CLSIM_DATACLASSES_ONLY YES)
endif (OPENCL_FOUND)

if (APPLE AND (CMAKE_SYSTEM_VERSION VERSION_LESS "11.0.0"))
  colormsg(RED   "+-- OpenCL is not compatible with clsim in this version of OS X. Only building the clsim dataclasses.")
  set (BUILD_CLSIM_DATACLASSES_ONLY YES)
  set (tmp_OPENCL_FOUND ${OPENCL_FOUND})
  set(OPENCL_FOUND False)
endif (APPLE AND (CMAKE_SYSTEM_VERSION VERSION_LESS "11.0.0"))


# The IceTray project definition
if (BUILD_CLSIM_DATACLASSES_ONLY)
  i3_project(clsim
    DOCS_DIR resources/docs
  )

  if (OPENCL_FOUND)
    colormsg(RED   "+-- Only building the clsim dataclasses. (manual override, OpenCL *is* available)")
  else (OPENCL_FOUND)
    colormsg(RED   "+-- OpenCL is not installed on your system. Only building the clsim dataclasses.")
  endif (OPENCL_FOUND)
else (BUILD_CLSIM_DATACLASSES_ONLY)
  i3_project(clsim
    PYTHON_DIR python
    DOCS_DIR resources/docs
  )
endif (BUILD_CLSIM_DATACLASSES_ONLY)


SET(LIB_${PROJECT_NAME}_SOURCEFILES
  # private/clsim
  private/clsim/I3CLSimEventStatistics.cxx
  private/clsim/I3Photon.cxx
  private/clsim/I3CompressedPhoton.cxx
  private/clsim/I3CLSimFlasherPulse.cxx
  private/clsim/util/I3MuonSlicer.cxx
  private/clsim/util/I3MuonSliceRemoverAndPulseRelabeler.cxx
  private/clsim/util/I3TauSanitizer.cxx
  private/clsim/dom/I3PhotonToMCPEConverter.cxx
  private/clsim/shadow/I3ShadowedPhotonRemover.cxx
  private/clsim/shadow/I3ShadowedPhotonRemoverModule.cxx
  private/clsim/shadow/I3ExtraGeometryItem.cxx
  private/clsim/shadow/I3ExtraGeometryItemUnion.cxx
  private/clsim/shadow/I3ExtraGeometryItemMove.cxx
  private/clsim/shadow/I3ExtraGeometryItemCylinder.cxx

  # tableio converters
  private/clsim/converter/I3PhotonConverter.cxx
  private/clsim/converter/I3MCHitConverterWithIDs.cxx
  private/clsim/converter/I3MCPEConverterWithIDs.cxx
)

SET(LIB_${PROJECT_NAME}_PROJECTS icetray dataio dataclasses phys-services tableio simclasses)
SET(LIB_${PROJECT_NAME}_TOOLS python boost gsl)


# check for Geant4
if (GEANT4_FOUND AND OPENCL_FOUND AND NOT DISABLE_GEANT4_IN_CLSIM)
    ADD_DEFINITIONS(-DHAS_GEANT4)

    LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
        private/geant4/TrkCerenkov.cxx
        private/geant4/TrkDetectorConstruction.cxx
        private/geant4/TrkEMPhysicsUHE.cxx
        private/geant4/TrkEventAction.cxx
        private/geant4/TrkOpticalPhysics.cxx
        private/geant4/TrkPrimaryGeneratorAction.cxx
        private/geant4/TrkStackingAction.cxx
        private/geant4/TrkUISessionToQueue.cxx
        private/geant4/TrkUserEventInformation.cxx
    )
    LIST(APPEND LIB_${PROJECT_NAME}_TOOLS clhep geant4)

    colormsg(GREEN "+-- Geant4 support")
elseif (GEANT4_FOUND AND OPENCL_FOUND AND DISABLE_GEANT4_IN_CLSIM)
    colormsg(YELLOW "+-- Geant4 is installed on your system but has been explicitly disabled using -DDISABLE_GEANT4_IN_CLSIM. clsim will fail if it is not used with parameterizations.")
else (GEANT4_FOUND AND OPENCL_FOUND AND NOT DISABLE_GEANT4_IN_CLSIM)
    colormsg(RED   "+-- Geant4 is not installed on your system. clsim will fail if it is not used with parameterizations.")
endif (GEANT4_FOUND AND OPENCL_FOUND AND NOT DISABLE_GEANT4_IN_CLSIM)


# check for numpy. The C++ part of the tabulator needs it.
if(NUMPY_FOUND)
  ADD_DEFINITIONS(-DUSE_NUMPY)

  LIST(APPEND LIB_${PROJECT_NAME}_TOOLS numpy)
  ADD_DEFINITIONS(-DUSE_TABULATOR)
  SET(BUILD_TABULATOR YES)
  colormsg(GREEN "+-- numpy support (for tabulator)")
else(NUMPY_FOUND)
  SET(BUILD_TABULATOR NO)
  colormsg(GREEN "+-- no numpy support (tabulator will not be built)")
endif(NUMPY_FOUND)

# The icecube master header file I3.h seems to include inttypes.h without setting
# __STDC_FORMAT_MACROS. clsim makes use of the PRI* preprocessor definitions,
# so it needs those macros. Add the definition right here to work around this.
ADD_DEFINITIONS(-D__STDC_FORMAT_MACROS)


# check for Q-Frame support
if(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)
  # we assume that the existence of this module means that we have
  # to deal with a Q-frame-enabled meta-project
  ADD_DEFINITIONS(-DIS_Q_FRAME_ENABLED)
  colormsg(GREEN "+-- Q-frame-enabled")
else(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)
  colormsg(YELLOW "+-- old-style, non-Q-frame module")
endif(EXISTS ${CMAKE_SOURCE_DIR}/dataio/private/dataio/QConverter.cxx)


# check for portable binary archive in "icetray" (instead of cmake/boost)
if(EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/portable_binary_archive.hpp)
  ADD_DEFINITIONS(-DHAS_PBA_IN_ICETRAY)
  colormsg(GREEN "+-- portable_binary_archive is in the icetray project")
else(EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/portable_binary_archive.hpp)
  colormsg(YELLOW "+-- portable_binary_archive is *not* in the icetray project")
endif(EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/portable_binary_archive.hpp)


# Check for the existence of icetray/python/dataclass_suite.hpp
# If this file does not exist, a local copy of the python
# binding helper headers is used. TODO: This needs to go away
# as soon as an offline-software release contains these headers.
# (They already are in icetray/trunk.)
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/python/dataclass_suite.hpp)
  ADD_DEFINITIONS(-DNO_PYTHON_DATACLASS_SUITE)
  colormsg(YELLOW "+-- icetray does not have dataclass_suite.hpp. Using local version.")
endif(NOT EXISTS ${CMAKE_SOURCE_DIR}/icetray/public/icetray/python/dataclass_suite.hpp)


# This is the core clsim stuff (everything using OpenCL and/or Geant4)
if(NOT BUILD_CLSIM_DATACLASSES_ONLY)

  LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES
    # private/clsim
    private/clsim/I3CLSimMediumProperties.cxx
    private/clsim/I3CLSimModule.cxx
    private/clsim/I3CLSimModuleHelper.cxx
    private/clsim/I3CLSimLightSourceParameterization.cxx
    private/clsim/I3CLSimLightSourceToStepConverterPPC.cxx
    private/clsim/I3CLSimLightSourceToStepConverter.cxx
    private/clsim/I3CLSimLightSourceToStepConverterUtils.cxx
    private/clsim/I3CLSimPhoton.cxx
    private/clsim/I3CLSimPhotonHistory.cxx
    private/clsim/random_value/I3CLSimRandomValueApplyFunction.cxx
    private/clsim/random_value/I3CLSimRandomValue.cxx
    private/clsim/random_value/I3CLSimRandomValueHenyeyGreenstein.cxx
    private/clsim/random_value/I3CLSimRandomValueInterpolatedDistribution.cxx
    private/clsim/random_value/I3CLSimRandomValueMixed.cxx
    private/clsim/random_value/I3CLSimRandomValueRayleighScatteringCosAngle.cxx
    private/clsim/random_value/I3CLSimRandomValueSimplifiedLiu.cxx
    private/clsim/random_value/I3CLSimRandomValueWlenCherenkovNoDispersion.cxx
    private/clsim/random_value/I3CLSimRandomValueNormalDistribution.cxx
    private/clsim/random_value/I3CLSimRandomValueFixParameter.cxx
    private/clsim/random_value/I3CLSimRandomValueConstant.cxx
    private/clsim/random_value/I3CLSimRandomValueUniform.cxx
    private/clsim/I3CLSimSimpleGeometryFromI3Geometry.cxx
    private/clsim/I3CLSimSimpleGeometryTextFile.cxx
    private/clsim/I3CLSimSimpleGeometryUserConfigurable.cxx
    private/clsim/I3CLSimStep.cxx
    private/clsim/function/I3CLSimFunctionAbsLenIceCube.cxx
    private/clsim/function/I3CLSimFunctionConstant.cxx
    private/clsim/function/I3CLSimFunctionDeltaPeak.cxx
    private/clsim/function/I3CLSimFunction.cxx
    private/clsim/function/I3CLSimFunctionFromTable.cxx
    private/clsim/function/I3CLSimFunctionPolynomial.cxx
    private/clsim/function/I3CLSimFunctionRefIndexIceCube.cxx
    private/clsim/function/I3CLSimFunctionRefIndexQuanFry.cxx
    private/clsim/function/I3CLSimFunctionScatLenIceCube.cxx
    private/clsim/function/I3CLSimFunctionScatLenPartic.cxx
    private/clsim/function/I3CLSimScalarField.cxx
    private/clsim/function/I3CLSimScalarFieldConstant.cxx
    private/clsim/function/I3CLSimScalarFieldAnisotropyAbsLenScaling.cxx
    private/clsim/function/I3CLSimScalarFieldIceTiltZShift.cxx
    private/clsim/function/I3CLSimVectorTransform.cxx
    private/clsim/function/I3CLSimVectorTransformConstant.cxx
    private/clsim/function/I3CLSimVectorTransformMatrix.cxx
    private/clsim/I3CLSimLightSource.cxx
    private/clsim/I3CLSimSpectrumTable.cxx
    private/clsim/I3CLSimLightSourceToStepConverterFlasher.cxx

    # private/geant4
    private/geant4/I3CLSimLightSourceToStepConverterGeant4.cxx

    # private/opencl/
    private/opencl/I3CLSimHelperGenerateGeometrySource.cxx
    private/opencl/I3CLSimHelperGenerateMediumPropertiesSource.cxx
    private/opencl/I3CLSimHelperGenerateMediumPropertiesSource_Optimizers.cxx
    private/opencl/I3CLSimStepToPhotonConverterOpenCL.cxx
    private/opencl/I3CLSimOpenCLDevice.cxx
    private/opencl/ieeehalfprecision.cxx

    # private/test/
    private/test/I3CLSimMediumPropertiesTester.cxx
    private/test/I3CLSimRandomDistributionTester.cxx
    private/test/I3CLSimTesterBase.cxx
    private/test/I3CLSimFunctionTester.cxx
    private/test/I3CLSimScalarFieldTester.cxx
    private/test/I3CLSimVectorTransformTester.cxx
  )

  #LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS )
  LIST(APPEND LIB_${PROJECT_NAME}_TOOLS opencl)

  if(EXISTS $ENV{I3_DATA}/safeprimes_base32.txt)
    colormsg(CYAN   "+-- $ENV{I3_DATA}/safeprimes_base32.txt data file exists, skipping download")
  elseif(NOT EXISTS ${CMAKE_SOURCE_DIR}/clsim/resources/safeprimes_base32.txt)
    colormsg(RED   "+-- Downloading and unpacking the safeprimes_base32.txt data file...")
    colormsg(RED   "    (this file will only be downloaded once)")
    EXECUTE_PROCESS(COMMAND ${CMAKE_SOURCE_DIR}/clsim/resources/download_safeprimes.sh ${CMAKE_SOURCE_DIR}/clsim/resources)
    colormsg(GREEN "+-- ...download complete.")
  else(NOT EXISTS ${CMAKE_SOURCE_DIR}/clsim/resources/safeprimes_base32.txt)
    colormsg(CYAN  "+-- safeprimes_base32.txt already downloaded")
  endif(EXISTS $ENV{I3_DATA}/safeprimes_base32.txt)
endif(NOT BUILD_CLSIM_DATACLASSES_ONLY)

# the clsim library definition
i3_add_library(${PROJECT_NAME}
  ${LIB_${PROJECT_NAME}_SOURCEFILES}
  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
  USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
  )

# the make-safeprimes tool needs gmp, so only compile it if that tool is available
if(GMP_FOUND)
  i3_executable(make_safeprimes
    private/make_safeprimes/main.cxx
    USE_TOOLS python boost gmp
    )
  i3_executable(convert_safeprimes
    private/make_safeprimes/convert.cxx
    USE_PROJECTS icetray
    USE_TOOLS python boost
    )
  colormsg(GREEN "+-- gmp support (make_safeprimes utility)")
else(GMP_FOUND)
  colormsg(CYAN  "+-- no gmp support (make_safeprimes utility)")
endif(GMP_FOUND)


add_subdirectory(private/pybindings)

## restore OPENCL_FOUND
set(OPENCL_FOUND ${tmp_OPENCL_FOUND})
