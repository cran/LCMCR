\name{lcmCR}
\alias{lcmCR}
\title{
Bayesian Nonparametric Latent Class Capture-Recapture
}
\description{
Create and initialize an object of class \code{lcm_CR_Basic}.
}
\usage{
lcmCR(captures, tabular = FALSE, in_list_label = "1", not_in_list_label = "0", 
  K = 5, a_alpha = 0.25, b_alpha = 0.25, buffer_size = 10000, thinning = 10, 
  seed = "auto", verbose = TRUE)
}
\arguments{
  \item{captures}{
	input dataset. A data frame with the multiple-recapture data. See 'Details' for input formats.
  }
  \item{tabular}{
	a logical value indicating whether or not the data is tabulated. See 'Details'.
  }
  \item{in_list_label}{
	factor label that indicates that individual is in list (e.g. 'Yes')
  }
  \item{not_in_list_label}{
	factor label that indicates that individual is in not list (e.g. 'No')
  }
  \item{K}{
	maximum number of latent classes. Indicates the truncation level of the stick-breaking process.
  }
  \item{a_alpha}{
	shape parameter of  the prior distribution of concentration parameter of the stick-breaking process.
  }
  \item{b_alpha}{
	inverse scale parameter of  the prior distribution of concentration parameter of the stick-breaking process.
  }
  \item{buffer_size}{
	size of the tracing buffer.
  }
  \item{thinning}{
	thinning interval for the tracing buffer
  }
  \item{seed}{
	integer seed of the internal RNG. 
  }
  \item{verbose}{
	Generate progress messages?
  }
}
\details{
	Input data must be provided as a data frame.  The first J columns are two-level factors representing the multiple-recapture lists.  Arguments \code{in_list_label} and \code{not_in_list_label} indicate the labels that represent inclusion and exclusion from the lists.  This function supports two input formats:
\itemize{
	\item When \code{tabular=FALSE} each row represents a single individual's capture history. The number of rows must match  the size of the observed population. Rows indicating no capture in all list simultaneously are illegal.
	\item When \code{tabular=TRUE} each row represents a particular capture pattern. This format requires an additional numeric column at the right indicating the count corresponding to such pattern. 
}
}
\value{
An object of class \code{\link{lcm_CR_Basic}} initialized, after running \code{burnin} iterations.
}
%\references{
%}
\author{
Daniel Manrique-Vallier
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lcm_CR_Basic}}, \code{\link{lcm_CR_Basic_generator}}
}
\examples{
require('LCMCR')
data(kosovo_aggregate)
sampler <- lcmCR(captures = kosovo_aggregate, tabular = FALSE, in_list_label = '1',
  not_in_list_label = '0', K = 10, a_alpha = 0.25, b_alpha = 0.25,
  seed = 'auto', buffer_size = 10000, thinning = 100)
sampler
N <- lcmCR_PostSampl(sampler, burnin = 10000, samples = 1000, thinning = 100)
quantile(N, c(0.025, 0.5, 0.975))
}

\keyword{ capture-recapture }
