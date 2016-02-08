\name{LCMCR-package}
\alias{LCMCR-package}
\alias{LCMCR}
\docType{package}
\title{
Capture-Recapture Estimation using Bayesian Nonparametric latent-class models}
\description{
This package implements a fully Bayesian multiple-recapture method for estimating the unknown size of a population using non-parametric latent class models. This is an implementation of the method described in Manrique-Vallier (forthcoming). The estimation algorithm is based on Markov Chain Monte Carlo sampling. 
}
\details{
\tabular{ll}{
Package: \tab LCMCR\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2016-02-07\cr
License: \tab GPL >= 2\cr
}
 }
\author{
Daniel Manrique-Vallier \email{dmanriqu@indiana.edu}
}
\references{
Manrique-Vallier, D. (forthcoming) "Bayesian Population Size Estimation Using Dirichlet Process Mixtures", Biometrics. 
}

\examples{
library('LCMCR')

###Using Kosovo data.###
data(kosovo_aggregate)

###Example 1: Capture-Recapture estimation using convenience functions###
#Create and initialize an LCMCR object for MCMC sampling#
sampler <- lcmCR(captures = kosovo_aggregate, tabular = FALSE, in_list_label = '1',
      not_in_list_label = '0', K = 10, a_alpha = 0.25, b_alpha = 0.25, 
      seed = 'auto', buffer_size = 10000, thinning = 100)
#Obtain 1000 samples from the posterior distribution of N#
N <- lcmCR_PostSampl(sampler, burnin = 10000, samples = 1000, thinning = 100)

#Posterior quantiles#
quantile(N, c(0.025, 0.5, 0.975))

###Example 2: Capture-Recapture estimation using the lcm_CR_Basic object directly###
#Create and initialize an LCMCR object for MCMC sampling#
sampler <- lcmCR(captures = kosovo_aggregate, tabular = FALSE, in_list_label = '1',
                 not_in_list_label = '0', K = 10, a_alpha = 0.25, b_alpha = 0.25, 
                 seed = 'auto', buffer_size = 1000, thinning = 100)

#Run 10000 iterations as burn-in
sampler$Update(10000)

#List all parameters from the model
sampler$Get_Param_List()

#Set parameter 'n0' for tracing
sampler$Set_Trace('n0')

#List currently traced parameters.
sampler$Get_Trace_List()

#Activate tracing
sampler$Activate_Tracing()

#Run the sampler 100000 times
sampler$Update(100000)

#Get the 1000 samples from the posterior distribution of N
N <- sampler$Get_Trace('n0') + sampler$n

#Plot the trace of N
plot(N, type = 'l')

#Compute posterior quantiles
quantile(N, c(0.025, 0.5, 0.975))

}
