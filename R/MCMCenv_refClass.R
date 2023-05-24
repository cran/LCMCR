 # 
 # Copyright (C) 2012-2023 Daniel Manrique-Vallier
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

##########################################
# Functions and objects definitions for 
#  manipulationg MCMC environments.
# (implementation with reference classes)

#this is a base class. Has to be used by inheritance. NEED TO SET THE POINTER!!!
#To do: VALIDATION!

setRefClass(
  Class="MCMCenviron",
  fields=list(
    pointer='externalptr', 
    blobsize = 'numeric',
	local_seed ='numeric'
  ), 
  methods=list(
    initialize = function(){
      #.self$initialized = FALSE
		usingMethods('Init_Model', 
		'Update',
		'Get_Iteration',
		'Get_Trace_List', 
		'Get_Param_List',
		'Get_Trace',
		'Get_Param',
		'Get_Trace_Size',
		'Set_Trace',
		'Reset_Traces', 
		'Activate_Tracing', 
		'Deactivate_Tracing', 
		'Change_SubSamp', 
		'Change_Trace_Length', 
		'Set_Seed') 
    },
    Init_Model = function(output = TRUE, seed = c('auto','r.seed')){
	"Initializes the sampler. Set output = FALSE to suppress console output. Seed can be an integer, 'auto' or 'r.seed'."
	  #seed the internal RNG
	  if (!is.numeric(seed)){
			seed <- match.arg(seed)
	  }
	  if (seed == 'r.seed'){
		runif(1)
		.self$local_seed <- Reduce(bitwXor, .Random.seed)
		.self$Set_Seed(.self$local_seed)
	  } else if (is.numeric(seed)){
		.self$local_seed <- as.integer(seed)
		.self$Set_Seed(.self$local_seed)
	  } #default behavior: do nothing. C++ class will seed automatically from clock.
	  #Control messages
	  if (!interactive()) output = FALSE
	  if (output){
		tmp <- .Call('R_Activate_Chain_Messages', .self$pointer)
	  } else {
		tmp <- .Call('R_Deactivate_Chain_Messages', .self$pointer)
	  }
      #SEXP R_Init_Model(SEXP p);
      tmp <- .Call('R_Init_Model', .self$pointer)
    },
    Update = function(iter, output = TRUE){
	"Runs num_iter iterations of the sampler. Set output = FALSE to suppress console output."
      #SEXP R_Update_Model(SEXP p, SEXP int_iter)
	  if (!interactive()) output = FALSE
	  if (output){
		tmp <- .Call('R_Activate_Updating_Output', .self$pointer)
	  } else {
		tmp <- .Call('R_Deactivate_Updating_Output', .self$pointer)
	  }
      tmp <- .Call('R_Update_Model', .self$pointer, as.integer(iter))
    },
    Get_Iteration = function(){
	"Retrieves the current number of iterations the sampler."
      return(.Call('R_Get_Iteration', .self$pointer))
    },
    Get_Trace_List = function(){
	"Retrieves the names of the parameters being currently traced."
      tmp <- .Call('R_Get_Trace_List', .self$pointer)
      return(tmp)
    },
    Get_Param_List = function(){
	"Retrieves the names of the parameters of the model."
      tmp <- .Call('R_Get_Param_List', .self$pointer)
      return(tmp)
    },
    Get_Trace = function(param){
	"Retrieves the tracing buffer for parameter param. Returns an array whose first dimension is the sample index."
      tmp <- .Call('R_Get_Trace', .self$pointer, as.character(param))
      return(aperm(tmp, length(dim(tmp)):1))
    },
    Get_Param = function(param){
	"Retrieves the current value of the parameter param."
      tmp <- .Call('R_Get_Param', .self$pointer, as.character(param))
      return(aperm(tmp, length(dim(tmp)):1))
    },
   Get_Trace_Size = function(){
   "Retrieves the size (in iterations) of the trace buffer."
      tmp <- .Call('R_Get_Trace_Size', .self$pointer)
      return(tmp)
    },
    Set_Trace = function(params){
	"Adds parameters to the tracer. To list the available parameters for tracing use the Get_Param_List() method."
      for (tra in params){
        tmp <- .Call('R_Set_Trace', .self$pointer, as.character(tra))
      }
    },
    Reset_Traces = function(){
	"Deletes the content of the tracing buffer."
      tmp <- .Call('R_Reset_Traces', .self$pointer)
    },
    Activate_Tracing = function(){
	"Activates the tracing buffer."
      tmp <- .Call('R_Activate_Tracing', .self$pointer)
    },
    Deactivate_Tracing = function(){
	"Deactivates the tracing buffer."
      tmp <- .Call('R_Deactivate_Tracing', .self$pointer)
    },
    Change_SubSamp = function(new_subsamp = 1){
	"Changes the sub-sampling period (thinning) of the tracing buffer."
      tmp <- .Call('R_Change_SubSamp', .self$pointer, as.integer(new_subsamp))
    },
    Change_Trace_Length = function(new_length){
	"Changes the size (in number of samples) of the tracing buffer."
      tmp <- .Call('R_Change_Trace_Size', .self$pointer, as.integer(new_length))
    },
	Set_Seed = function(seed){
	"Seeds the internal random number generator. Doesn't affect R's internal RNG."
		#SEXP R_Set_Seed(SEXP p, SEXP seed);
		tmp <- .Call('R_Set_Seed', .self$pointer, as.integer(seed))
		.self$local_seed <- seed
	},
	Get_Status = function(){
	#SEXP R_Get_Status(SEXP p)
	#Returns a vector c(<iteration>, <initialized>, <buffer size>, <buffer used>, <tracer activated>, <thinning>)  		
		tmp <- .Call('R_Get_Status', .self$pointer)
		return(list(
			iteration = tmp[1],
			initialized = (tmp[2] == 1),
			buffer_size = tmp[3],
			buffer_used = tmp[4],
			tracing = (tmp[5] == 1),
			thinning = tmp[6]
			)
		)	
	},
	show = function(){
		tmp <- .self$Get_Status()
		cat('\tInitialized =', tmp$initialized,'\n')
		cat('\tCurrent iteration =', tmp$iteration,'\n')
		cat('Tracer:\n')
		cat('\tActivated =', tmp$tracing, '\n')
		cat('\tCapacity =',tmp$buffer_size,'samples.\n')
		cat('\tUsed  = ', tmp$buffer_used, ' (', round(tmp$buffer_used/tmp$buffer_size*100, digits = 2), '%)\n', sep ='');
		cat('\tThining =', tmp$thinning,'\n')
		cat('\tCurrently Tracing:', paste(.self$Get_Trace_List(), collapse=', '), '\n');
	})
  )


