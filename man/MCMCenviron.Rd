\name{MCMCenviron-class}
\Rdversion{1.1}
\docType{class}
\alias{MCMCenviron-class}

\title{Class \code{"MCMCenviron"}}
\description{
A generic interface for MCMC sampler objects implementing Bayesian models. Methods provide access to underlying functionality implemented in C++. The underlying implementation provides basic functionality for controlling the chain, and a 'tracing buffer' for storing and retrieving the samples.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}

\author{
Daniel Manrique-Vallier
}
\note{
This class is not designed to be used directly, but as a generic interface for samplers implementing specific models.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("MCMCenviron")
}
\section{Fields}{
(All fields are read-only.)
  \describe{
    \item{\code{pointer}:}{external pointer to the C++ object }
	\item{\code{blobsize}:}{size (in bytes) of the raw object data for serialization. (currently not implemented.)}
    \item{\code{seed}:}{seed of the internal random number generator. }
  }
}
\section{Methods}{
GENERAL METHODS
  \describe{
    \item{\code{Init_Model(output = TRUE, seed=c('auto', 'r.seed'))}:}{ 
	Initializes the sampler.
	\subsection{Arguments}{
			\describe{
				\item{\code{output:}}{logical. Print messages to the screen?}
				\item{\code{seed:}}{integer. Seed of the internal RNG. Additionally, \code{seed='auto'} autogenerates the seed from the internal clock; \code{seed='r.seed'} autogenerates the seed from the current state of the \code{.Random.seed} variable.}
			}
		}
	}
    \item{\code{Update(num_iter, output = TRUE)}:}{ Runs \code{num_iter} iterations of the sampler. Set \code{output = FALSE} to suppress console output.}
    \item{\code{Get_Iteration()}:}{ Retrieves the current number of iterations the sampler. }
    \item{\code{Get_Param_List()}:}{ Retrieves the names of the parameters of the model.}
    \item{\code{Get_Param(param)}:}{ Retrieves the current value of the parameter \code{param}. }
    \item{\code{Set_Seed(seed)}:}{ Seeds the internal random number generator. It does not affect R's internal RNG. }
	\item{\code{Get_Status():}}{
		Retrieves the current state of the chain
		\subsection{Value}{
			\describe{
				\item{\code{iteration}}{numeric. Current iteration} 
				\item{\code{initialized}}{logical. Is the sampler initialized?}
				\item{\code{buffer_size}}{numeric. Capacity (in samples) of the tracing buffer.}
				\item{\code{buffer_used}}{numeric. Number of samples currently stored in the tracing buffer.}
				\item{\code{tracing}}{character. Names of the variables currently traced.}
				\item{\code{thinning}}{numeric. Thinning interval of the tracing buffer.}
			}
		}
	}
  }
  
  METHODS FOR CONTROLLING THE TRACING BUFFER
  \describe{
    \item{\code{Get_Trace_List()}:}{ Retrieves the names of the parameters being currently traced.}
    \item{\code{Activate_Tracing()}:}{ Activates the tracing buffer. Traced variables will be stored in the buffer when generated with \code{Update()}.}
    \item{\code{Deactivate_Tracing()}:}{ Deactivates the tracing buffer. Calls to \code{Update()} will not store samples in the buffer.}
    \item{\code{Set_Trace(traces)}:}{
		Adds parameters to tracer. 
		\subsection{Arguments}{
			\describe{
				\item{\code{param}: }{ character vector. Names of the parameters to trace. To list the available parameters for tracing use the \code{Get_Param_List()} method.}
			}
		}
	}
    \item{\code{Get_Trace(param)}:}{ Retrieves samples stored in the tracing buffer.
		\subsection{Arguments}{
			\describe{
				\item{\code{param}: }{ character. Name of the parameter to retrieve.}
			}
		}
		\subsection{Value}{An array. The first dimension indexes the sample; the rest correspond to the original dimensions of the parameter as defined in the model.
		}
	}
    \item{\code{Reset_Traces()}:}{ Deletes the content of the tracing buffer. }
    \item{\code{Change_SubSamp(new_subsamp)}:}{ Changes the sub-sampling period (thinning) of the tracing buffer. 
		\subsection{Warning}{
			This operation deletes the current content of the tracing buffer.
		}
	}
    \item{\code{Get_Trace_Size()}:}{ Retrieves the size (in number of samples) of the trace buffer. }
    \item{\code{Change_Trace_Length(new_length)}:}{ Changes the size (in number of samples) of the tracing buffer. 
		\subsection{Warning}{
			This operation deletes the current content of the tracing buffer.
		}
	}
}
}
