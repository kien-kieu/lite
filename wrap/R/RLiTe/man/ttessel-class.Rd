% Line Tessellation (LiTe) library
% Release 1.1
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% Copyright INRA 2006-2015.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{TTessel}
\alias{TTessel}
\title{C++ T-tessellation class}
\description{C++ class exposed to R through Rcpp. A TTessel object
  represents a planar T-tessellation of a rectangular domain.}
\section{Documentation}{
  Methods exposed to R can be listed using \code{show}. Documentation of
  input arguments can be found in LiTe Doxygen documentation.}
\examples{
## load a LineTes object from file
ltes <- new(LineTes)
ltes$read(file.path(path.package("RLiTe"),"tessellation_data.txt"))
## generate a TTessel object from the LineTes object
myTes <- new(TTessel,ltes)
plot(myTes,asp=1)
## other available methods
show(TTessel)
myTes$getCellAreas()
myTes$getInternalSegmentNumber()
myTes$number_of_blocking_segments()
}