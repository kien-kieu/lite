% Line Tessellation (LiTe) library
% Release 1.1
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% Copyright INRA 2006-2015.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{LineTes}
\alias{LineTes}
\title{C++ LineTes class}
\description{C++ class exposed to R through Rcpp. A LineTes object
  represents a planar polygonal tessellation of a rectangular domain.}
\section{Documentation}{
  Methods exposed to R can be listed using \code{show}. Documentation of
  input arguments can be found in LiTe Doxygen documentation.}
\examples{
show(LineTes)
tes <- new(LineTes)
tes$read(file.path(path.package("RLiTe"),"tessellation_data.txt"))
plot(tes,asp=1)
}