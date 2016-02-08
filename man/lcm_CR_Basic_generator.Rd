\name{lcm_CR_Basic_generator}
\alias{lcm_CR_Basic_generator}

\title{Generator for Class \code{lcm_CR_Basic}}
\description{
Generator function for class \code{\linkS4class{lcm_CR_Basic}}. 
}
\usage{
lcm_CR_Basic_generator(...)
}
\arguments{
 \item{...}{ arguments to be passed to \code{\linkS4class{lcm_CR_Basic}} constructor.}
}
	
\value{An object of class  \code{\linkS4class{lcm_CR_Basic}}.}
\author{
Daniel Manrique-Vallier.}
\note{
The convenience function \code{\link{lcmCR}} provides a simpler mechanism to create \code{lcm_CR_Basic} objects.
}

\seealso{
\code{\link{lcmCR}}
}
\examples{
data(kosovo_aggregate)
x <- lcm_CR_Basic_generator(data_captures=kosovo_aggregate, K=10, a_alpha=0.25, b_alpha=0.25, 
                                len_buffer=10000, subsamp=500, in_list_symbol = '1')
x$Get_Status()
}
