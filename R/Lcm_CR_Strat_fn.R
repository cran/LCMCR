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
# R5 class definition for lcm_CR_Strat 
# Depends on functions in: 
# - ../R_General/ArrayUtils.R
# - ../R_General/MCMCenv_refClass.j
# - ../R_Interface/CR_Support.R
# (c) Daniel Manrique-Vallier 2019
######################################################################

lcm_CR_Strat_generator <- setRefClass(
  Class = "lcm_CR_Strat",
  fields = list( J = 'numeric', K = 'numeric', C = 'numeric', n = 'numeric', 
                 nC = 'table', Captures = 'data.frame', Stratification = 'factor' ),
  methods = list(
    initialize = function(
      strata, data_captures, K, 
      a_alpha = 0.25, b_alpha = 0.25, a_lambda = 4, b_lambda = 0.1, 
      min_n = 4, min_lists = 2, N_Max_factor = 5,
      in_list_symbol = '1', len_buffer, subsamp)
    {
      callSuper()
      zeros <- apply(as.matrix(data_captures != in_list_symbol), MARGIN = 1, FUN = prod)==1
      if (sum(zeros) > 0){
        data_captures <- data_captures[!zeros,]
        strata <- strata[!zeros]
        warning('Inconsistent rows. ', sum(zeros), ' "no capture" rows eliminated.')
      }
      .self$J = NCOL(data_captures)
      .self$n = NROW(data_captures)
      .self$C = nlevels(strata)
      .self$K = K
      .self$Captures = data_captures
      .self$Stratification = addNA(strata)
      .self$nC = table (.self$Stratification)
      #SEXP R_Create_LCM_CR_Strat(
      #	SEXP x_flat, SEXP n_strat_plus_NA, SEXP J_plus_strat,
      #	SEXP n, SEXP K,
      #	SEXP a_alpha, SEXP b_alpha, SEXP a_lambda, SEXP b_lambda, 
      #	SEXP min_n, SEXP min_lists, SEXP N_Max_factor,
      #	SEXP len_buffer, SEXP subsamp
      #);
      tmp <- .Call('R_Create_LCM_CR_Strat', 
                   as.integer(fn_factor2CR(df.factors = .self$Captures, 
                                           sym.capture = in_list_symbol,
                                           strat = .self$Stratification )
                              ), 
                   as.integer(.self$C + 1), # Number of strata + NA
                   as.integer(.self$J + 1), # J + stratification
                   as.integer(.self$n), 
                   as.integer(.self$K), 
                   as.double(a_alpha), as.double(b_alpha), 
                   as.double(a_lambda), as.double(b_lambda), 
                   as.integer(min_n), as.integer(min_lists), as.double(N_Max_factor),
                   as.integer(len_buffer), as.integer(subsamp)
      )
      .self$pointer = tmp
    }
  ),
  contains = "MCMCenviron"
)


lcmCR_strat <- function(captures, tabular = FALSE, in_list_label = '1', not_in_list_label = '0', K = 5, a_alpha = 0.25, b_alpha=0.25, buffer_size=10000, thinning = 10, seed = 'auto', verbose = TRUE){
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

lcmCR_Strat_PostSampl <- function(object, burnin=10000, samples = 1000, thinning = 10, clear_buffer = FALSE, output = TRUE){
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
