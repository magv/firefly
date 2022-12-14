cmake_minimum_required(VERSION 3.1)
project(FireFly VERSION 2.0.3 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# decleare some colours
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(FLINT REQUIRED)

# set the build type to release. If one needs debugging symbols, compile
# with "-DCMAKE_BUILD_TYPE=Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable all warnings in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# set paths
set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source)
add_definitions("-DFIREFLY_ROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
set(INCLUDE_PATH ${SOURCE_PATH}/include)
set(INCLUDE_PATH_SRC ${SOURCE_PATH}/include/firefly)

# include headers
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories(${INCLUDE_PATH})
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${FLINT_INCLUDE_DIRS})
set(PKG_INCLUDE -I${FLINT_INCLUDE_DIRS} -I${ZLIB_INCLUDE_DIRS})

set(MPI_FLAG "0")

set(ENABLE_SHARED TRUE CACHE BOOL "Build and install the shared library")
set(ENABLE_STATIC TRUE CACHE BOOL "Build and install the static library")
set(ENABLE_FF_INSERT TRUE CACHE BOOL "Build and install the 'ff_insert' executable")
set(ENABLE_FIREFLY FALSE CACHE BOOL "Build and install the 'firefly' executable")
set(ENABLE_EXAMPLE TRUE CACHE BOOL "Build the 'example' executable")
set(WITH_MPI FALSE CACHE BOOL "Build with MPI support")


# Check for MPI
if(${WITH_MPI})
  set(MPI_FLAG "1")
  if("${MPI_CXX_INCLUDE_PATH}" STREQUAL "")
    find_package(MPI REQUIRED)
  else()
    message("-- Using user-defined MPI implementation library: ${MPI_CXX_LIBRARIES}")
    message("-- Using user-defined MPI implementation headers: ${MPI_CXX_INCLUDE_PATH}")
  endif()
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(PKG_INCLUDE ${PKG_INCLUDE} -I${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif()

# Firefly library files
set(LIBSOURCES
  ${SOURCE_PATH}/AmplitudeParser.cpp
  ${SOURCE_PATH}/BaseReconst.cpp
  ${SOURCE_PATH}/FFInt.cpp
  ${SOURCE_PATH}/FFIntVec.cpp
  ${SOURCE_PATH}/FFThieleInterpolator.cpp
  ${SOURCE_PATH}/gzstream.cpp
  ${SOURCE_PATH}/HornerGenerator.cpp
  ${SOURCE_PATH}/Monomial.cpp
  ${SOURCE_PATH}/ParserUtils.cpp
  ${SOURCE_PATH}/Poly.cpp
  ${SOURCE_PATH}/Polynomial.cpp
  ${SOURCE_PATH}/PolynomialFF.cpp
  ${SOURCE_PATH}/PolyReconst.cpp
  ${SOURCE_PATH}/RationalFunction.cpp
  ${SOURCE_PATH}/RationalFunctionFF.cpp
  ${SOURCE_PATH}/RationalNumber.cpp
  ${SOURCE_PATH}/RatReconst.cpp
  ${SOURCE_PATH}/ReconstHelper.cpp
  ${SOURCE_PATH}/ShuntingYardParser.cpp
  ${SOURCE_PATH}/utils.cpp)

set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/firefly")
set(BIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INCLUDE_INSTALL_DIR_2 "${CMAKE_INSTALL_PREFIX}/include")
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")

configure_file (
  "${INCLUDE_PATH_SRC}/version.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/version.hpp")

configure_file (
  "${INCLUDE_PATH_SRC}/config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/config.hpp")

if(${ENABLE_STATIC})
  add_library(FireFly_static STATIC ${LIBSOURCES})
  set_target_properties(FireFly_static PROPERTIES OUTPUT_NAME firefly)
endif()
if(${ENABLE_SHARED})
  add_library(FireFly_shared SHARED ${LIBSOURCES})
  set_target_properties(FireFly_shared PROPERTIES OUTPUT_NAME firefly)
endif()

set(LIBS ${FLINT_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(${WITH_MPI})
  #string(STRIP ${MPI_CXX_LINK_FLAGS} MPI_CXX_LINK_FLAGS)
  #set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
  set(EXAMPLE_EXE ${PROJECT_SOURCE_DIR}/example_mpi.cpp)
else()
  set(EXAMPLE_EXE ${PROJECT_SOURCE_DIR}/example.cpp)
  set(INSERT_IBPS_EXE ${PROJECT_SOURCE_DIR}/ff_insert.cpp)
  set(FIREFLY_EXE ${PROJECT_SOURCE_DIR}/firefly.cpp)
endif()

if ("${WITH_JEMALLOC}" STREQUAL "true")
  execute_process(COMMAND "jemalloc-config" "--libdir" OUTPUT_VARIABLE JEMALLOC_LIB_DIR)
  execute_process(COMMAND "jemalloc-config" "--libs" OUTPUT_VARIABLE JEMALLOC_LIBS)
  string(REGEX REPLACE "\n$" "" JEMALLOC_LIB_DIR "${JEMALLOC_LIB_DIR}")
  string(REGEX REPLACE "\n$" "" JEMALLOC_LIBS "${JEMALLOC_LIBS}")
  if(NOT "${JEMALLOC_LIB_DIR}" STREQUAL "")
    message("-- Found jemalloc: ${JEMALLOC_LIB_DIR}")
    set(JEMALLOC_LIB "-L${JEMALLOC_LIB_DIR} -Wl,-rpath,${JEMALLOC_LIB_DIR} -ljemalloc" ${JEMALLOC_LIBS})
  else()
    message("-- jemalloc-config not found, searching for pkgconfig file instead")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JEMALLOC REQUIRED jemalloc)
    set(JEMALLOC_LIB "-L${JEMALLOC_LIBDIR} -Wl,-rpath,${JEMALLOC_LIBDIR} -ljemalloc")
  endif()
  set(LIBS ${JEMALLOC_LIB} ${LIBS})
endif()

if(${WITH_MPI})
  if(${ENABLE_EXAMPLE} AND ${ENABLE_STATIC})
    add_executable(example ${EXAMPLE_EXE})
    target_link_libraries(example FireFly_static ${LIBS})
  endif()
  if(${ENABLE_SHARED})
    target_link_libraries(FireFly_shared ${LIBS})
  endif()
else()
  if(${ENABLE_EXAMPLE} AND ${ENABLE_STATIC})
    add_executable(example ${EXAMPLE_EXE})
    target_link_libraries(example FireFly_static ${LIBS})
  endif()
  if(${ENABLE_SHARED})
    target_link_libraries(FireFly_shared ${LIBS})
  endif()
  if(${ENABLE_FF_INSERT} AND ${ENABLE_STATIC})
    add_executable(ff_insert ${INSERT_IBPS_EXE})
    target_link_libraries(ff_insert FireFly_static ${LIBS})
  endif()
  if(${ENABLE_FIREFLY} AND ${ENABLE_STATIC})
    message("-- Building executable `firefly'. Note that it is still in development and not officially released.")
    add_executable(firefly ${FIREFLY_EXE})
    target_link_libraries(firefly FireFly_static ${LIBS})
  endif()
endif()

if("${DISABLE_INSTALL}" STREQUAL "")
  # install libraries
  if(${ENABLE_STATIC})
    install(
      TARGETS FireFly_static
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib)
  endif()
  if(${ENABLE_SHARED})
    install(
      TARGETS FireFly_shared
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib)
  endif()

  set(PKG_LIBS ${LIBS})
  list(REMOVE_DUPLICATES PKG_LIBS)
  string(REPLACE ";" " " PKG_LIBS "${PKG_LIBS}")

  set(PKG_INCLUDE ${PKG_INCLUDE})
  list(REMOVE_DUPLICATES PKG_INCLUDE)
  string(REPLACE ";" " " PKG_INCLUDE "${PKG_INCLUDE}")

  configure_file("cmake/pc.in" "firefly.pc" @ONLY)

  # install pkgconfig file
  install(FILES
    "${CMAKE_BINARY_DIR}/firefly.pc"
    DESTINATION "${LIB_INSTALL_DIR}/pkgconfig")

  # install public headers
  install(FILES
    "${INCLUDE_PATH_SRC}/AmplitudeParser.hpp"
    "${INCLUDE_PATH_SRC}/BaseReconst.hpp"
    "${INCLUDE_PATH_SRC}/BlackBoxBase.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/config.hpp"
    "${INCLUDE_PATH_SRC}/DenseSolver.hpp"
    "${INCLUDE_PATH_SRC}/FFInt.hpp"
    "${INCLUDE_PATH_SRC}/FFIntVec.hpp"
    "${INCLUDE_PATH_SRC}/FFThieleInterpolator.hpp"
    "${INCLUDE_PATH_SRC}/gzstream.hpp"
    "${INCLUDE_PATH_SRC}/HornerGenerator.hpp"
    "${INCLUDE_PATH_SRC}/Logger.hpp"
    "${INCLUDE_PATH_SRC}/Monomial.hpp"
    "${INCLUDE_PATH_SRC}/ParserUtils.hpp"
    "${INCLUDE_PATH_SRC}/Poly.hpp"
    "${INCLUDE_PATH_SRC}/Polynomial.hpp"
    "${INCLUDE_PATH_SRC}/PolynomialFF.hpp"
    "${INCLUDE_PATH_SRC}/PolyReconst.hpp"
    "${INCLUDE_PATH_SRC}/RationalFunction.hpp"
    "${INCLUDE_PATH_SRC}/RationalFunctionFF.hpp"
    "${INCLUDE_PATH_SRC}/RationalNumber.hpp"
    "${INCLUDE_PATH_SRC}/RatReconst.hpp"
    "${INCLUDE_PATH_SRC}/ReconstHelper.hpp"
    "${INCLUDE_PATH_SRC}/Reconstructor.hpp"
    "${INCLUDE_PATH_SRC}/ShuntingYardParser.hpp"
    "${INCLUDE_PATH_SRC}/ThreadPool.hpp"
    "${INCLUDE_PATH_SRC}/tinydir.h"
    "${INCLUDE_PATH_SRC}/UintHasher.hpp"
    "${INCLUDE_PATH_SRC}/utils.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/include/firefly/version.hpp"
    DESTINATION "${INCLUDE_INSTALL_DIR}")

  # install ibp insertion executable
  if(NOT ${WITH_MPI})
    if(${ENABLE_FF_INSERT} AND ${ENABLE_STATIC})
      install(TARGETS ff_insert DESTINATION "${BIN_INSTALL_DIR}")
    endif()
    if(${ENABLE_FIREFLY} AND ${ENABLE_STATIC})
      install(TARGETS firefly DESTINATION "${BIN_INSTALL_DIR}")
    endif()
  endif()

  if(${WITH_MPI})
    install(FILES "${INCLUDE_PATH_SRC}/MPIWorker.hpp" DESTINATION "${INCLUDE_INSTALL_DIR}")
  endif()

  # adding documentation
  add_subdirectory(doc)

  # test
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
    enable_testing()
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
    add_subdirectory(test)
  endif()

  # benchmarks
  if("${BUILD_BENCH}" STREQUAL true)
    if(${ENABLE_STATIC})
      add_executable(benchmarks ${PROJECT_SOURCE_DIR}/benchmarks/benchmarks.cpp)
      target_link_libraries(benchmarks FireFly_static ${LIBS})
      add_executable(benchmarks_no_bt ${PROJECT_SOURCE_DIR}/benchmarks/benchmarks_no_bt.cpp)
      target_link_libraries(benchmarks_no_bt FireFly_static ${LIBS})
    endif()
  endif()

  if(${ENABLE_SHARED})
    add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue --bold --no-newline "FireFly info: ")
    add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E echo "To run the example type: ./example")
    add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue --bold --no-newline "FireFly info: ")
    add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E echo "To run the inserter type: ./ff_insert")
    add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue --bold --no-newline "FireFly info: ")
    add_custom_command(TARGET FireFly_shared POST_BUILD COMMAND ${CMAKE_COMMAND} -E echo "To compile and run tests type: make check")
  endif()
endif()
