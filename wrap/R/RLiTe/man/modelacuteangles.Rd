% Line Tessellation (LiTe) library
% Release 1.1
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% Copyright INRA 2006-2015.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{ModelAcuteAngles}
\alias{ModelAcuteAngles}
\title{Energy for Gibbsian T-tessellations with penalty on angles
  between segments}
\description{Create an Energy object representing a Gibbsian
  T-tessellation with penalty on the acute angles between segments}
\usage{ModelAcuteAngles(tau,beta)}
\arguments{
  \item{tau}{chemical activity parameter}
  \item{beta}{parameter controling the strength of penalization}
}
\details{
  The energy has the following form
  \deqn{
    -\stackrel{\circ}{n}_{\mathrm{s}}(T)\log\tau
    + \beta\sum_{v}(\frac{\pi}{2}-\phi(v))
  }{
    -n(T)*log(tau) + beta*sum(pi/2-phi(v))
  }
  where
  \itemize{
    \item \eqn{\stackrel{\circ}{n}_{\mathrm{s}}(T)}{n(T)} is the number of
    internal segments of the T-tessellation \eqn{T},
    \item \eqn{\sum_{v}}{sum()} is a sum taken over all tessellation
    vertices,
    \item \eqn{\phi(v)}{phi(v)} is the acute angle between the two
    segments meeting at vertex \eqn{v}.
  }
}
\value{An Energy object}
\examples{
myTes <- new(TTessel)
myTes$setDomain(1,1)
myMod <- ModelAcuteAngles(tau=12.1,beta=2.5)
myMod$set_ttessel(myTes)
}

