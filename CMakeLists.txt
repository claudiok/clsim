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


i3_add_library(clsim
  private/clsim/*.cxx
  private/geant4/*.cxx
  private/opencl/*.cxx
  USE_TOOLS python boost log4cplus clhep geant4 opencl
  USE_PROJECTS icetray dataio dataclasses phys-services
  )

i3_executable(make_safeprimes
  private/make_safeprimes/main.c
  #USE_PROJECTS g4sim
  USE_TOOLS python boost gmp
  )

add_subdirectory(private/pybindings)

