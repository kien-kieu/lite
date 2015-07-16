# Line Tessellation (LiTe) library
# Release 1.1
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# Copyright INRA 2006-2015.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

AC_INIT([src/rlite_module.cpp])
## Check for non-standard programs: lite-config(1)
AC_PATH_PROG([LITE_CONFIG], [lite-config])
## If lite-config was found, let's use it
if test "${LITE_CONFIG}" != ""; then
   # Use lite-config for linker arguments
   LITE_PREFIX=`${LITE_CONFIG} --prefix`
   LITE_LIBDIR=${LITE_PREFIX}/lib
   LITE_INCLUDEDIR=${LITE_PREFIX}/include
else
   AC_MSG_ERROR([lite-config not found, is LiTe installed?])
fi
# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(LITE_LIBDIR) 
AC_SUBST(LITE_INCLUDEDIR) 
AC_OUTPUT(src/Makevars)

