\name{lcmCR_PostSampl}
\alias{lcmCR_PostSampl}
\title{
Generate Samples from the Posterior Distribution of Population Size under a LCMCR Model
}
\description{
Convenience function for generate samples from the posterior distribution of the population size using an initialized \code{\link{lcm_CR_Basic}} object.
}
\usage{
lcmCR_PostSampl(object, burnin = 10000, samples = 1000, thinning = 10, 
    clear_buffer = FALSE, output = TRUE)
}
\arguments{
\item{object}{
an initialized \code{\link{lcm_CR_Basic}} object.
}
\item{burnin}{
number of burn in iterations.
}
\item{samples}{
Nnmber of samples to be generated.  Note that this is not the same as the number of iterations for the sampler.  Samples are saved one every \code{thinning} iterations.
}
\item{thinning}{
subsampling interval.  Samples are saved one every \code{thinning} iterations.
}
\item{clear_buffer}{
logical. Clear the tracing buffer before sampling?
}
\item{output}{
logical. Print messages?
}
}
\value{
A vector with the \code{samples} posterior samples of the population size parameter.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Daniel Manrique-Vallier
}
\note{
To create and initialize the lcm_CR_Basic object use \code{\link{lcmCR}} or \code{\link{lcm_CR_Basic_generator}}.  The user is responsible to check whether the chain has reached the stationary distribution or not.
}

\section{Warning}{
Invoking this function deletes the content of the object's tracing buffer.
}


\examples{
data(kosovo_aggregate)
sampler <- lcmCR(captures = kosovo_aggregate, tabular = FALSE, in_list_label = '1',
      not_in_list_label = '0', K = 10, a_alpha = 0.25, b_alpha = 0.25, seed = 'auto')
N <- lcmCR_PostSampl(sampler, burnin = 10000, samples = 1000, thinning = 100)
quantile(N, c(0.025, 0.5, 0.975))
}

