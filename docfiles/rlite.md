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

Download the package source at https://github.com/kien-kieu/lite/blob/release-x.y/build/wrap/R/RLiTe_x.y.tar.gzwhere x.y is the (latest) release number of LiTe.


Installation of RLiTe
---------------------
RLiTe can be installed as a standard R package. For instance, the following command line should do the job.

    R CMD INSTALL RLiTe_x.y.tar.gz

Or within a R interactive session, RLiTe can be installed running

   install.packages("RLiTe_x.y.tar.gz",repos=NULL)

Windows users
-------------
Both CGAL and Rcpp are available under Windows. So what's the problem? To make them talk! The only successful way of installing RLiTe under Windows we found until now is based on [MSYS2](http://sourceforge.net/projects/msys2). The procedure below was set up by Rémy Drouilhet.

First install MSYS2. A standard install is OK.

Edit the PATH environment variable, add msysroot\mingw32\bin to PATH where msysroot is the directory where MSYS2 has been installed.

Next open the MSYS2 shell and install the CGAL package

    pacman -S mingw-w64-i686-cgal 

Finally download the Windows binaries versions of Rcpp and RLite available from LiTe download pages and install them (first Rcpp, then RLiTe).

The above procedure may not be very convenient if you are using Rcpp for other R packages. We are trying to regularly update the Rcpp binaries provided on LiTe download pages, but it may happen to be slightly outdated. In such a case, you may build yourself the Windows binaries for Rcpp. Download the Rcpp package source and build the binary following the R documentation guidelines.