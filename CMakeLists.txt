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

include("FindOpenCL.cmake")

if(OPENCL_FOUND)
  message(STATUS "  + OpenCL includes at ${OPENCL_INCLUDE_DIRS}")
  message(STATUS "  + OpenCL libraries at ${OPENCL_LIB_DIR}")
  message(STATUS "  + OpenCL library at ${OPENCL_LIBRARIES}")
  if( OPENCL_HAS_CPP_BINDINGS )
      message(STATUS "  + OpenCL has CPP bindings. Full include is: " ${OPENCL_INCLUDE_DIRS} )
  else( OPENCL_HAS_CPP_BINDINGS )
      message(STATUS "  + No OpenCL CPP bindings found" )
  endif( OPENCL_HAS_CPP_BINDINGS )

  include_directories(${OPENCL_INCLUDE_DIRS})
else(OPENCL_FOUND)
  message(STATUS "  - OpenCL *not* found")
  set(OPENCL_LIBRARIES "") 
  add_definitions(-DI3CLSIM_WITHOUT_OPENCL)
endif(OPENCL_FOUND)



i3_add_library(clsim
  private/clsim/*.cxx
  private/geant4/*.cxx
  private/opencl/*.cxx
  USE_TOOLS python boost log4cplus clhep geant4
  USE_PROJECTS icetray dataio dataclasses phys-services

  LINK_LIBRARIES ${OPENCL_LIBRARIES}
  )


## this is a hack and should not be here. move to cmake project
#set(GMP_VERSION "4.3.2")
#tooldef(gmp
#  include/gmp-${GMP_VERSION}
#  gmp/gmp.h
#  lib/gmp-${GMP_VERSION}
#  NONE
#  gmp
#  )
#i3_executable(make_safeprimes
#  private/make_safeprimes/main.c
#  USE_TOOLS python boost gmp
#  )

add_subdirectory(private/pybindings)

