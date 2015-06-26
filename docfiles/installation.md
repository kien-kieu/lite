Installation
============

<!-- Line Tessellation (LiTe) library
     |||Development version
     Authors: Katarzyna Adamczyk and Kiên Kiêu.
     |||Copyright INRA 2006-yyyy.
     Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
     License: GPL v3. -->

A prerequisite for installing LiTe is the installation of CGAL (http://www.cgal.org) and CMake 2.8.x. The generation of documentation requires Doxygen (http://www.doxygen.org).

Download
--------

### Release version of LiTe
Source files can be downloaded as an archive from https://mulcyber.toulouse.inra.fr/frs/?group_id=xx. Unpack the archive into a given directory.

### Latest version of LiTe
Source files of the latest development version of LiTe are also available through git:

     git clone http://mulcyber.toulouse.inra.fr/anonscm/git/lite/lite.git

Build LiTe
----------

The installation of LiTe is based on CMake. There are several ways to build LiTe using either the command line interface or graphical CMake clients. The few command lines below are provided as a simple example of building. The path name `lite_src` refers to the directory where the LiTe source files have been unpacked.

    cd lite_src
    mkdir ../lite_build
    cd ../lite_build
    cmake ../lite_src

There are a few CMake options specific to LiTe:

- BUILD_DOCUMENTATION. OFF (default) or ON. If ON, Doxygen documentation is generated. Of course, Doxygen must be installed. 
- STANDALONE_RLITE. OFF (default) or ON. If ON, a standalone version of the R package RLiTe is generated which does not depend on the LiTe library.

As an example, if documentation should be generated, build LiTe using the command line below.

    cmake ../lite_src -D BUILD_DOCUMENTATION=ON

Compile and install LiTe
------------------------
This stage depends on your development environment. For Linux users, it is performed using the standard commands

    make
    make install

Libraries are copied into CMAKE_INSTALL_PREFIX/lib, header files into CMAKE_INSTALL_PREFIX/include, HTML documentation files into CMAKE_INSTALL_PREFIX/share/doc/lite. A shell script lite-config is copied to CMAKE_INSTALL_PREFIX. Compiled demos and applications are left in the build directory. 

Of course, installation requires write access to CMAKE_INSTALL_PREFIX. If needed, run

    sudo make install

The shell script lite-config can be used to get the installation prefix a posteriori.

Compile of a new application using LiTe
---------------------------------------

Compiling programs using LiTe may require additional steps depending on the OS. For instance, under Linux, if CMAKE_INSTALL_PREFIX is /usr/local (default value), one may need to run ldconfig as root. Or if CMAKE_INSTALL_PREFIX is set to some non standard path say "foo", one should have the environment variable LD_LIBRARY_PATH containing "foo/lib".

Compilation can be done using the template CMake configuration file `CMakeLists_standalone_template.txt` provided in the build directory.

Uninstalling LiTe
-----------------

There is no user-friendly way of uninstalling LiTe. Under a Unix-like OS and if the build directory is still there, installed filed can be removed using the command line

    xargs rm < install_manifest.txt
