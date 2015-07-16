# Line Tessellation (LiTe) library
# Release 1.1
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# Copyright INRA 2006-2015.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

project( liteApplications )
cmake_minimum_required(VERSION 2.8.0)
find_package(CGAL QUIET COMPONENTS Core )
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  add_executable(checkACS checkACS)
  target_link_libraries(checkACS ttesselStaticLib)
# count STL has not been updated after major changes on arrangement classes in
# CGAL. It cannot be compiled any more.
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

