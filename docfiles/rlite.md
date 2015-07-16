The R package RLiTe {#rlite}
===================

<!-- Line Tessellation (LiTe) library
     |||Development version
     Authors: Katarzyna Adamczyk and Kiên Kiêu.
     |||Copyright INRA 2006-yyyy.
     Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
     License: GPL v3. -->

LiTe main features are also available under [R](http://www.r-project.org) as a package named RLiTe. RLiTe is based on LiTe and the Rcpp package (we became aware of the latter thanks to Rémy Drouilhet).

Prerequisites
-------------
The following must be installed:
- CGAL library (http://www.cgal.org).
- R package Rcpp.

Download
--------
At the moment only the installation from package source is really supported. This means in practice that mostly Unix-like systems (Linux, Mac OS X) are supported. No standard Windows binaries are provided. The procedure for Windows given at the end of this pages is still rather experimental.

Download the package source at https://github.com/kien-kieu/lite/blob/release-x.y/build/wrap/R/RLiTe_x.y.tar.gz where x.y is the (latest) release number of LiTe.


Installation of RLiTe
---------------------
RLiTe can be installed as a standard R package. For instance, the following command line should do the job.

    R CMD INSTALL RLiTe_x.y.tar.gz

Or within a R interactive session, RLiTe can be installed running

   install.packages("RLiTe_x.y.tar.gz",repos=NULL)

Windows users
-------------
Both CGAL and Rcpp are available under Windows. So what's the problem? To make them talk! The only successful way of installing RLiTe under Windows we found until now is based on [MSYS2](http://msys2.github.io). The procedure below was set up by Rémy Drouilhet.

First install MSYS2. A standard install is OK. You may install either the 32-bits or the 64-bits version (assumed below). *However, at the moment, it is only possible to install the 32-bits version of RLiTe.* Follow the instructions provided on MSYS2 Web pages in order to achieve the installation of MSYS2.

Next open the MinGW-w64 Win32 terminal and install gcc, make and CGAL:

     pacman -S mingw-w64-i686-gcc
     pacman -S make
     pacman -S mingw-w64-i686-cgal

Download the RLiTe package source (see above).

Reopen the Mingw-w64 Win32 terminal, run R (32-bits version), e.g.:

       /c/R/R-3.2.1/bin/i386/R --no-save

and install both Rcpp and RLiTe

       install.packages("Rcpp",repos="http://cran.univ-lyon1.fr",INSTALL_opts="--no-multiarch")
       install.packages("RLiTe_x.y.tar.gz",INSTALL_opts="--no-multiarch")
       
Add to the PATH environment variable of Windows the absolute path to the subdirectory mingw32\bin of MSYS2. For instance, the path may be C:\\msys64\\mingw32\\bin.

Run R as usual (32-bits version) and load RLiTe.