% Line Tessellation (LiTe) library
% Release 1.1
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% Copyright INRA 2006-2015.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{ModelAreasAngles}
\alias{ModelAreasAngles}
\title{Energy for Gibbsian T-tessellations with penalty on squared areas
and acute angles}
\description{Create an Energy object representing a Gibbsian
  T-tessellation with penalty on the sum of squared cell areas and of
  angles between segments}
\usage{ModelAreasAngles(tau,alpha,beta)}
\arguments{
  \item{tau}{chemical activity parameter}
  \item{alpha}{parameter controling the strength of the penalization on
    squared cell areas}
  \item{beta}{parameter controling the strength of the penalization on
    angles between segments}
}
\details{
  The energy has the following form
  \deqn{
    -\stackrel{\circ}{n}_{\mathrm{s}}(T)\log\tau
    + \alpha\sum_{c}a(c)^2 + \beta\sum_{v}(\frac{\pi}{2}-\phi(v))
  }{
    -n(T)*log(tau) + alpha*sum(a(cell)^2) + beta*sum(pi/2-phi(v))
  }
  Notation is explained in documentation of functions
  \code{\link{ModelSquaredAreas}} and \code{\link{ModelAcuteAngles}}.
}
\value{An Energy object}
\examples{
myTes <- new(TTessel)
myTes$setDomain(1,1)
myMod <- ModelAreasAngles(tau=2,alpha=93000,beta=200)
myMod$set_ttessel(myTes)
}

