% Line Tessellation (LiTe) library
% Release 1.1
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% Copyright INRA 2006-2015.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{getSMFSamplingPeriod}
\alias{getSMFSamplingPeriod}
\title{Compute a sampling period of an SMF chain based on segment
  renewal}
\description{Support for sampling an SMF chain (simulation of Gibbsian
  T-tessellations. The proposed sampling period is determined based on
  the renewal rate of tessellation segments.}
\usage{getSMFSamplingPeriod(sim,mod,tes,target.renew)}
\arguments{
  \item{sim}{An SMF simulator as an Rcpp_SMFChain object}
  \item{mod}{A Gibbsian model of T-tessellation as an Rcpp_Energy
    object}
  \item{tes}{The T-tessellation handled by sim as an Rcpp_TTessel
    object}
  \item{target.renew}{The target renewal rate of tessellation
  segments. Must be between 0 and 1.}
}
\value{The proposed sampling period as a number of iterations}
\details{
  A T-tessellation is made of segments (maximal sets of contiguous and
  aligned edges). When simulating a Gibbsian T-tessellation, segments
  are born, change their size and die. Given a time interval, there are
  segments living at the beginning of the interval. A fraction of them
  die before the end of the interval. The renewal rate of segments is
  defined as the latter fraction. The aim of the getSMFSamplingPeriod is
  to determine how long must be the time interval in order to achieve a
  given renewal rate.

  The method implemented in getSMFSamplingPeriod is rather empirical and
  approximative. The key idea is to run the simulator (SMF chain) during
  some time and to register the segment ages. In particular the output
  is random.
}
\examples{
myTes <- new(TTessel)
myTes$setDomain(1,1)
myMod <- ModelCRTT(tau=1.9)
myMod$set_ttessel(myTes)
mySim <- new(SMFChain,myMod,0.33,0.33)
setLiTeSeed(5)
mySim$step(500) 
getSMFSamplingPeriod(mySim,myMod,myTes,0.75)
}
