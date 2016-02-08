\name{kosovo_aggregate}
\alias{kosovo_aggregate}
\docType{data}
\title{
Killings in the Kosovo war from March 20 to June 22, 1999.
}
\description{
Capture pattern data for $J = 4$ independently collected lists that jointly document $n = 4400$ observed killings in the Kosovo war betwen March 20 to June 22, 1999.
}
\usage{data("kosovo_aggregate")}
\format{
  A data frame with 4400 observations on the following 4 variables.
  \describe{
    \item{\code{EXH}}{a factor with levels \code{0} \code{1}}
    \item{\code{ABA}}{a factor with levels \code{0} \code{1}}
    \item{\code{OSCE}}{a factor with levels \code{0} \code{1}}
    \item{\code{HRW}}{a factor with levels \code{0} \code{1}}
  }
}
\details{
This data set was analyzed by Ball et al. (2002).}
%\source{
%
%}
\references{
Ball, P., Betts, W., Scheuren, F., Dudukovic, J., and Asher, J. (2002), ``Killings and
Refugee Flow in Kosovo, MarchJune, 1999," Report to ICTY.
}
\examples{
data(kosovo_aggregate)
}
\keyword{datasets}
