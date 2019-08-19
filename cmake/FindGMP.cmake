# Copyright (c) 2006, Laurent Montel, <montel@kde.org>
# Copyright (c) 2007, Francesco Biscani, <bluescarni@gmail.com>
# Copyright (c) 2019, Jonas Klappert and Fabian Lange

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ------------------------------------------------------------------------------------------

# Try to find the GMP libraries used with C++ code:
# GMP_FOUND - System has GMP lib
# GMP_INCLUDE_DIR - The GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP

#if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
  # Force search at every time, in case configuration changes
#  unset(GMP_INCLUDE_DIR CACHE)
#  unset(GMP_LIBRARIES CACHE)
#endif (GMP_INCLUDE_DIR AND GMP_LIBRARIES)

find_path(GMP_INCLUDE_DIR NAMES gmpxx.h)

set(GMP_FIND_VERSION_MAJOR 6)
set(GMP_FIND_VERSION_MINOR 1)
set(GMP_FIND_VERSION_PATCH 2)
set(GMP_FIND_VERSION
    "${GMP_FIND_VERSION_MAJOR}.${GMP_FIND_VERSION_MINOR}.${GMP_FIND_VERSION_PATCH}")

if(GMP_INCLUDE_DIR)
  # Since the GMP version macros may be in a file included by gmp.h of the form
  # gmp-.*[_]?.*.h (e.g., gmp-x86_64.h), we search each of them.
  file(GLOB GMP_HEADERS "${GMP_INCLUDE_DIR}/gmp.h" "${GMP_INCLUDE_DIR}/gmp-*.h")
  foreach(gmp_header_filename ${GMP_HEADERS})
    file(READ "${gmp_header_filename}" _gmp_version_header)
    string(REGEX MATCH
      "define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" _gmp_major_version_match
      "${_gmp_version_header}")
    if(_gmp_major_version_match)
      set(GMP_MAJOR_VERSION "${CMAKE_MATCH_1}")
      string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)"
        _gmp_minor_version_match "${_gmp_version_header}")
      set(GMP_MINOR_VERSION "${CMAKE_MATCH_1}")
      string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
        _gmp_patchlevel_version_match "${_gmp_version_header}")
      set(GMP_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")
      set(GMP_VERSION
        ${GMP_MAJOR_VERSION}.${GMP_MINOR_VERSION}.${GMP_PATCHLEVEL_VERSION})
    endif()
  endforeach()

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT GMP_VERSION)
    set(GMP_VERSION_OK FALSE)
    message(STATUS "GMP version was not detected")
  elseif(${GMP_VERSION} VERSION_LESS ${GMP_FIND_VERSION})
    set(GMP_VERSION_OK FALSE)
  else()
    set(GMP_VERSION_OK TRUE)
  endif()
endif()

if(STBIN)
  find_library(GMP_LIBRARIES NAMES libgmp.a gmp)
else(STBIN)
  find_library(GMP_LIBRARIES NAMES libgmp.so gmp)
endif(STBIN)

if(GMP_INCLUDE_DIR AND GMP_LIBRARIES)
  set(GMP_FOUND TRUE)
endif(GMP_INCLUDE_DIR AND GMP_LIBRARIES)

if(GMP_FOUND AND GMP_VERSION_OK)
  message(STATUS "Found GMP ${GMP_VERSION} headers: ${GMP_INCLUDE_DIR}")
  message(STATUS "GMP library: ${GMP_LIBRARIES}")
else()
  message(FATAL_ERROR "GMP version ${GMP_VERSION} found in ${GMP_INCLUDE_DIR}, "
                   "but at least version ${GMP_FIND_VERSION} is required")
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)
