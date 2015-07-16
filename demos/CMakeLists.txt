# Line Tessellation (LiTe) library
# Release 1.1
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# Copyright INRA 2006-2015.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

project( liteDemos )
cmake_minimum_required(VERSION 2.8.0)
find_package(CGAL QUIET COMPONENTS Core )
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  add_executable(simACS sim_acs)
  # add_executable(plsa pseudo_lik_inference)
  # add_executable(exact pseudo_lik_exact)
  # add_executable(plareas pseudo_lik_areas)
  # add_executable(pldareas pseudo_lik_discr_areas)
  # add_executable(check_pldareas check_pseudo_discrete_value)
  # add_executable(check_pseudo check_pseudo)
  target_link_libraries(simACS ttesselStaticLib)
  # target_link_libraries(plsa ttesselStaticLib)
  # target_link_libraries(exact ttesselStaticLib)
  # target_link_libraries(plareas ttesselStaticLib)
  # target_link_libraries(pldareas ttesselStaticLib)
  # target_link_libraries(check_pldareas ttesselStaticLib)
  # target_link_libraries(check_pseudo ttesselStaticLib)
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

