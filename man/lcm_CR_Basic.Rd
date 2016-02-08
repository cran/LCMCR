\name{lcm_CR_Basic-class}
\Rdversion{1.1}
\docType{class}
\alias{lcm_CR_Basic-class}
\alias{lcm_CR_Basic}

\title{Class \code{"lcm_CR_Basic"}}
\description{
MCMC sampler for the Bayesian non-parametric latent class capture-recapture model.
}
\section{Extends}{
Class \code{"\linkS4class{MCMCenviron}"}, directly.  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
Daniel Manrique-Vallier
}
\note{
Use the convenience function \code{\link{lcmCR}} to create objects of this class.  This class inherits most of its functionality from \code{"\linkS4class{MCMCenviron}"}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lcmCR}}, \code{\linkS4class{MCMCenviron}}.
}
\examples{
showClass("lcm_CR_Basic")
}
\keyword{classes}
\section{Fields}{
	All fields are read-only.
  \describe{
	\item{\code{pointer}:}{external pointer to the C++ object.}
	\item{\code{blobsize}:}{size (in bytes) of the raw object data for serialization. (currently not implemented.)}
    \item{\code{local_seed}:}{seed of the internal random number generator. }
    \item{\code{J}:}{number of lists in the Capture-Recapture data.}
    \item{\code{K}:}{maximum number of latent classes in the model (truncation level of the stick-breaking process).}
    \item{\code{n}:}{observed number of individuals.}
    \item{\code{Captures}:}{original provided data.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{initialize(data_captures, K, a_alpha, b_alpha, in_list_symbol, len_buffer, subsamp)}:}{ 
	Class constructor.
		\subsection{Arguments}{
			\describe{
				  \item{\code{data_captures:}}{
					input dataset. A data frame with the multiple-recapture data. 
				  }
				  \item{\code{K}:}{
					maximum number of latent classes. Indicates the truncation level of the stick-breaking process.
				  }
				  \item{\code{a_alpha}:}{
					shape parameter of  the prior distribution of concentration parameter of the stick-breaking process.
				  }
				  \item{\code{b_alpha}:}{
					inverse scale parameter of  the prior distribution of concentration parameter of the stick-breaking process.
				  }
				  \item{\code{in_list_symbol}:}{
					factor label that indicates that individual is in list (e.g. 'Yes')
				  }
				  \item{\code{buffer_size}:}{
					Size of the tracing buffer.
				  }
				  \item{\code{subsamp}:}{
					thinning interval for the tracing buffer.
				  }
				  \item{\code{verbose}:}{
					logical. Generate progress messages?
				  }
			}
		}
	}
  }

The following methods are inherited (from the corresponding class):
Change_SubSamp ("MCMCenviron"), Set_Trace ("MCMCenviron"), Change_Trace_Length ("MCMCenviron"), initialize ("MCMCenviron"), Get_Iteration ("MCMCenviron"), Get_Param ("MCMCenviron"), Reset_Traces ("MCMCenviron"), Get_Status ("MCMCenviron"), Update ("MCMCenviron"), Get_Trace_Size ("MCMCenviron"), Get_Trace ("MCMCenviron"), Get_Trace_List ("MCMCenviron"), Get_Param_List ("MCMCenviron"), Init_Model ("MCMCenviron"), Activate_Tracing ("MCMCenviron"), Deactivate_Tracing ("MCMCenviron"), Set_Seed ("MCMCenviron"), show ("MCMCenviron")
}
