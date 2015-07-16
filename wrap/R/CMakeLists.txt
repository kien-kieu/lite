# Line Tessellation (LiTe) library
# Release 1.1
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# Copyright INRA 2006-2015.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

# Generation of the R package RLiTe
find_program(R_EXEC R)
if(R_EXEC)
  option(STANDALONE_RLITE "Generate a standalone version of RLiTe not depending on LiTe" OFF)
  add_custom_target(rlite_files COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/RLiTe ${CMAKE_CURRENT_BINARY_DIR}/RLiTe)
  # target lite_files is used for a standalone version of RLiTe not
  # depending on the installation of the LiTe library. Target building
  # consists in copying the source files of LiTe in the src directory
  # of RLiTe. All files in RLiTe/src related to autoconf are removed,
  # together with the script cleanup.
  add_custom_target(lite_files DEPENDS rlite_files COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/../../include/ttessel.h ${CMAKE_CURRENT_BINARY_DIR}/RLiTe/src/ttessel.h COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/ttessel.cpp ${CMAKE_CURRENT_BINARY_DIR}/RLiTe/src/ttessel.cpp COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/RLiTe/configure.ac COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/RLiTe/cleanup)
  string(TIMESTAMP TODAY "%Y-%m-%d")
  if(STANDALONE_RLITE)
    add_custom_target(rlite_prepare DEPENDS lite_files)
    set(PKGCXXFLAGS "")
    set(PKGLIBS "-lCGAL_Core -lCGAL -lgmpxx -lmpfr -lgmp -lboost_thread -lboost_system -lpthread")
    configure_file(rlite_makevars.in RLiTe/src/Makevars @ONLY)
    set(PKGLIBS_WIN "-lCGAL_Core -lCGAL -lgmpxx -lmpfr -lgmp -lboost_thread-mt -lboost_system-mt -lpthread")
    configure_file(RLiTe/src/Makevars.win.in RLiTe/src/Makevars.win @ONLY)
  else()
    set(PKGCXXFLAGS "-I@LITE_INCLUDEDIR@")
    set(PKGLIBS "-L@LITE_LIBDIR@ -lttesselSharedLib")
    add_custom_command(OUTPUT RLiTe/configure COMMAND aclocal
     COMMAND autoconf WORKING_DIRECTORY RLiTe)
    add_custom_target(configure_file DEPENDS rlite_files RLiTe/configure)
    add_custom_target(rlite_prepare DEPENDS configure_file)
    add_dependencies(rlite_prepare ttesselSharedLib)
    configure_file(rlite_makevars.in RLiTe/src/Makevars.in @ONLY)
  endif()
  configure_file(DESCRIPTION.in RLiTe/DESCRIPTION @ONLY)
  add_custom_target(rlite DEPENDS rlite_prepare COMMAND ${R_EXEC} CMD build RLiTe)
else()
  message(STATUS "Building of package RLiTe requires R")
endif()
