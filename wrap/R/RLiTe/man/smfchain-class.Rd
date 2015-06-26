% Line Tessellation (LiTe) library
% |||Development version
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% |||Copyright INRA 2006-yyyy.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{SMFChain}
\alias{SMFChain}
\title{C++ Gibbsian T-tessellation simulator}
\description{C++ class exposed to R through Rcpp. The Gibbsian
  T-tessellation model is simulated using a split, merge and flip Markov
chain.}
\section{C++ methods}{
  C++ methods are exposed through Rcpp. They are accessed using the
  \code{$} operator.
    \describe{
    \item{\code{setSMFProb(psplit,pmerge)}}{set the
      probabilities of proposing splits and merges. Flips are proposed
      with the complementary probability. Therefore the sum of both
      probabilities passed to \code{setSMFProb} must be less than 1.}
    \item{\code{step(n)}}{update the simulation Markov chain \code{n}
      times. Method \code{step} returns a matrix reporting the
      respective occurences of the three types of proposals (splits,
      merges and flips) and the numbers of accepted proposals.}
  }
}
\examples{
myTes <- new(TTessel)
myTes$setDomain(1,1)
myMod <- ModelCRTT(tau=1.9)
myMod$set_ttessel(myTes)
mySim <- new(SMFChain,myMod,0.33,0.33)
setLiTeSeed(5)
mySim$step(1000)
plot(myTes)
}