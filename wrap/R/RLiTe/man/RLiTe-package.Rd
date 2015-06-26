% Line Tessellation (LiTe) library
% |||Development version
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% |||Copyright INRA 2006-yyyy.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{RLiTe-package}
\alias{RLiTe-package}
\alias{RLiTe}
\docType{package}
\encoding{utf-8}
\title{
RLiTe is an embryonic interface to the C++ LiTe library
}
\description{
The purpose of LiTe is the simulation of Gibbs T-tessellations. RLiTe is an R interface based on Rcpp. At the moment, it is still at a very experimental stage.
}
\section{Main classes}{
  There are three main classes provided by RLiTe: TTessel, Energy and
  SMFChain. These classes are LiTe classes exposed to R via Rcpp. Use
  show method in order to list available methods:\cr\cr
  \code{show(TTessel)}
  }
\details{
\tabular{ll}{
Package: \tab RLiTe\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-27\cr
License: \tab LGPL\cr
}
}
\author{
Kiên Kiêu, Katarzyna Adamczyk

Maintainer: Kiên Kiêu <kien.kieu@jouy.inra.fr>
}
\references{
Kiêu, K., Adamczyk-Chauvat, K., Monod, H., Stoica. R. S. A completely random T-tessellation model and Gibbsian extensions. Spatial Statistics, 6, 118-138, 2013
}
\keyword{ package }
\keyword{ spatial }
\examples{
require(RLiTe)
myTes <- new(TTessel)
myTes$setDomain(2,1)
myMod <- ModelCRTT(tau=5.0)
myMod$set_ttessel(myTes)
mySim <- new(SMFChain,myMod,0.3,0.3)
setLiTeSeed(5)
mySim$step(100)
plot(myTes)
hist(myTes$getCellAreas())
}
