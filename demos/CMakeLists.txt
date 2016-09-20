# Line Tessellation (LiTe) library
# |||Development version
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# |||Copyright INRA 2006-yyyy.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

project( liteDemos )
cmake_minimum_required(VERSION 2.8.0)
find_package(CGAL QUIET COMPONENTS Core )
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  
  find_program(LITE_CONFIG lite-config)
  execute_process(COMMAND ${LITE_CONFIG} --prefix OUTPUT_VARIABLE
    LITE_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  include_directories(${LITE_PREFIX}/include)	
  include(${LITE_PREFIX}/lib/lite/lite.cmake)
  add_library(ttesselSharedLib SHARED IMPORTED)
  set_target_properties(ttesselSharedLib PROPERTIES IMPORTED_LOCATION
    ${LITE_PREFIX}/lib/libttesselSharedLib${CMAKE_SHARED_LIBRARY_SUFFIX})
  
  add_executable(simACS sim_acs)
  target_link_libraries(simACS ttesselStaticLib)

  add_executable(simACS2 sim_acs)
  target_link_libraries(simACS2 ttesselSharedLib)
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

