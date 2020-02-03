 # 
 # Copyright (C) 2012-2019 Daniel Manrique-Vallier
 # 
 # This program is free software; you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation; either version 2 of the License, or (at
 # your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful, but
 # WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 # General Public License for more details.
 # 
 # You should have received a copy of the GNU General Public License
 # along with this program; if not, write to the Free Software
 # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 #

######################################################################
# Functions for creating and manipulating NP_LCM objects
# Depends on functions in: 
# - ../R_General/ArrayUtils.R
# - ../R_General/MCMCenv_refClass.R
# - ../R_Interface/CR_Support.R
# (c) Daniel Manrique-Vallier 2019
######################################################################

#prototype:

lcm_CR_Basic_generator <- setRefClass(
  Class = "lcm_CR_Basic",
  fields = list(
	J = 'numeric',
	K = 'numeric',
	n = 'numeric',
	Captures = 'data.frame'
    ),
  methods = list(
    initialize =   function(data_captures, K, a_alpha, b_alpha, in_list_symbol = '1',
			len_buffer, subsamp){
      callSuper()
	  zeros <- apply(as.matrix(data_captures != in_list_symbol), MARGIN = 1, FUN = prod)==1
	  if (sum(zeros) > 0){
		data_captures <- data_captures[!zeros,]
		warning('Inconsistent rows. ', sum(zeros), ' "no capture" rows eliminated.')
	  }
      .self$J = NCOL(data_captures)
      .self$K = K
      .self$n = NROW(data_captures)
      .self$Captures = data_captures
	  #SEXP R_Create_LCM_CR_Basic(SEXP x_flat, SEXP J, SEXP n, SEXP K, SEXP Nmis_max, 
		#				SEXP a_alpha, SEXP b_alpha, 
		#				SEXP len_buffer, SEXP subsamp);
      tmp <- .Call('R_Create_LCM_CR_Basic', 
				   as.integer(fn_factor2CR(data_captures, in_list_symbol)), 
				   as.integer(.self$J), 
                   as.integer(.self$n), 
				   as.integer(K), 
                   as.double(a_alpha), 
				   as.double(b_alpha), 
                   as.integer(len_buffer), as.integer(subsamp)
      )
      .self$pointer = tmp
    }
  ),
  contains = "MCMCenviron"
)
.onUnload <- function (libpath) {
  gc(verbose = FALSE)
  library.dynam.unload("LCMCR", libpath)
}
lcmCR <- function(captures, tabular = FALSE, in_list_label = '1', not_in_list_label = '0', K = 5, a_alpha = 0.25, b_alpha=0.25, buffer_size=10000, thinning = 10, seed = 'auto', verbose = TRUE){
	if(is.matrix(captures)){
		captures <- fn_CRmatrix2dataframe(captures, in_list_label=in_list_label, not_in_list_label = not_in_list_label, tabulated = tabular)
	}
	if (tabular){
		captures <- fn_CRtabular2indiv(captures)
	}
	o <- lcm_CR_Basic_generator(data_captures = captures, K, a_alpha, b_alpha, in_list_symbol = in_list_label, len_buffer = buffer_size, subsamp = thinning)
	o$Init_Model(output=verbose, seed = seed)
	return(o)
}

lcmCR_PostSampl <- function(object, burnin=10000, samples = 1000, thinning = 10, clear_buffer = FALSE, output = TRUE){
  object$Update(burnin, output)
  object$Change_SubSamp(thinning)
  object$Change_Trace_Length(samples)
  if (!('n0' %in% object$Get_Trace_List())){
    object$Set_Trace('n0')
  }
  if (clear_buffer){
    object$Reset_Traces()
  }
  object$Activate_Tracing()
  object$Update(samples * thinning, output)
  N <- object$Get_Trace('n0') + object$n
  return(as.numeric(N))
}
