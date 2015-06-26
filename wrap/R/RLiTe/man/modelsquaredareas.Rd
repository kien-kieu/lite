% Line Tessellation (LiTe) library
% |||Development version
% Authors: Katarzyna Adamczyk and Kiên Kiêu.
% |||Copyright INRA 2006-yyyy.
% Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
% License: GPL v3.
\name{ModelSquaredAreas}
\alias{ModelSquaredAreas}
\title{Energy for Gibbsian T-tessellations with penalty on squared areas}
\description{Create an Energy object representing a Gibbsian
  T-tessellation with penalty on the sum of squared cell areas}
\usage{ModelSquaredAreas(tau,alpha)}
\arguments{
  \item{tau}{chemical activity parameter}
  \item{alpha}{parameter controling the strength of penalization}
}
\details{
  The energy has the following form
  \deqn{
    -\stackrel{\circ}{n}_{\mathrm{s}}(T)\log\tau
    + \alpha\sum_{c}a(c)^2
  }{
    -n(T)*log(tau) + alpha*sum(a(cell)^2)
  }
  where
  \itemize{
    \item \eqn{\stackrel{\circ}{n}_{\mathrm{s}}(T)}{n(T)} is the number of
    internal segments of the T-tessellation \eqn{T},
    \item \eqn{\sum_{c}a(c)^2}{sum(a(cell)^2)} is the sum of squared
    cell areas.
  }
}
\value{An Energy object}
\examples{
myTes <- new(TTessel)
myTes$setDomain(1,1)
myMod <- ModelSquaredAreas(tau=0.0043,alpha=10000)
myMod$set_ttessel(myTes)
}

