# Line Tessellation (LiTe) library
# |||Development version
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# |||Copyright INRA 2006-yyyy.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

project( lite )

cmake_minimum_required(VERSION 2.8.0)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  include_directories(BEFORE include)
  # Compilation of libraries, binaries and R package
  subdirs(include lib tests applications wrap/R)
  configure_file(lite-config.in lite-config @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/lite-config DESTINATION bin) 
  # Generation of the documentation
  option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
  if(BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    if (NOT DOXYGEN_FOUND)
      message(FATAL_ERROR
	"Doxygen is needed to build the documentation")
    endif()
    #-- Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in
      ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    #-- Add a custom target to run Doxygen when ever the project is built
    add_custom_target (Docs ALL
      COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
      SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
    install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION share/doc/lite
            PATTERN ".git" EXCLUDE)
  endif()
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()



