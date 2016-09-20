# Line Tessellation (LiTe) library
# |||Development version
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# |||Copyright INRA 2006-yyyy.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

project( lib_lite )
cmake_minimum_required(VERSION 2.8.9)
find_package(CGAL QUIET COMPONENTS Core )
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories (BEFORE ../include)

  add_library(ttesselSharedLib SHARED ttessel.cpp)
  # install(TARGETS ttesselSharedLib LIBRARY DESTINATION lib)
  install(TARGETS ttesselSharedLib DESTINATION lib)
  # property INSTALL_NAME_DIR is Mac OS X specific
  set_target_properties(ttesselSharedLib PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
  
  add_library(ttesselStaticLib STATIC ttessel)
  # set_target_properties(ttesselStaticLib PROPERTIES 
  #   CMAKE_POSITION_INDEPENDENT_CODE ON)
  install(TARGETS ttesselStaticLib  EXPORT lite ARCHIVE DESTINATION lib)
  install(EXPORT lite DESTINATION lib/lite)
 else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

