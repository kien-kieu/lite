% Line Tessellation (LiTe) library
% Release 1.1
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% Copyright INRA 2006-2015.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{ModelCRTT}
\alias{ModelCRTT}
\title{Energy for CRTT model}
\description{Create an Energy object representing a completely random
  T-tessellation}
\usage{ModelCRTT(tau=1.0)}
\arguments{
  \item{tau}{scaling parameter}
}
\details{
  The energy of the Completely Random T-Tessellation model has the form
  \deqn{
    -\stackrel{\circ}{n}_{\mathrm{s}}(T)\log\tau
  }{
    -n(T) log(tau)
  }
  where \eqn{\stackrel{\circ}{n}_{\mathrm{s}}(T)}{n(T)} is the number of
  internal segments of the T-tessellation \eqn{T}. The parameter
  \eqn{\tau}{tau} is a scaling parameter.
}
\value{An Energy object}
\examples{
myTes <- new(TTessel)
myTes$setDomain(1,1)
myMod <- ModelCRTT(tau=1.9)
myMod$set_ttessel(myTes)
}

