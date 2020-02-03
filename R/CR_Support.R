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
# (c) Daniel Manrique-Vallier 2019
######################################################################

.onUnload <- function (libpath) {
  gc(verbose = FALSE)
  library.dynam.unload("LCMCR", libpath)
}

fn_factor2CR <- function(df.factors, sym.capture = '1', missing = -1, offset = -1, 
                         C_style = TRUE, strat = NULL){
  #check that all variables have exactly 2 levels
  if (sum(fn_df_nlevels(df.factors) != 2) > 0 ) stop("Variable doesn't have 2 levels")
  if ( sum(sapply(fn_df_levels(df.factors), FUN = function(x) ! (sym.capture %in% x), simplify = T))  > 0)
    stop(paste( "'",sym.capture,"'", " not a level in all variables.", sep = ''  ))
  offset <- as.integer(offset); missing <- as.integer(missing)
  x <- ifelse(df.factors == sym.capture,as.integer(1),as.integer(0)) + offset + 1
  if (!is.null(strat)){
    s <- as.integer(as.integer(strat) + offset)
  } else {
    s <- NULL
  }
  r <- cbind(s, x) 
  r[is.na(r)] <- missing
  if(C_style){
    r <- t(r)
  }
  class(r) <- c('matrix', 'DMVmatrix')
  attr(r, 'levels') <- fn_df_levels(df.factors)
  attr(r, 'na.code') <- missing
  attr(r, 'offset') <- offset
  attr(r, 'C_style') <- C_style
  return(r)  
}
fn_CRfact2OnesZeros <- function(data_vector, in_list_label = 'Yes'){
  if(!is.factor(data_vector)) stop("Data vector is not a factor.")
  if(!(in_list_label %in% levels(data_vector))) stop(in_list_label," is not a level in data_vector")
  return(ifelse(data_vector == in_list_label, 1, 0))
}
fn_CRdataframe2matrix <- function(data, in_list_label = 'Yes'){
  if(!is.data.frame(data)) stop('Object "data" is not a dataframe')
  lastindx <- NCOL(data)
  res <- matrix(NA, ncol = NCOL(data), nrow = NROW(data))
  if(names(data)[lastindx] != 'Freq'){
    cols <- 1:lastindx
  } else {
    cols <- 1:(lastindx - 1)
    res[,lastindx] <- data[,lastindx] 
  }
  for (i in cols){
    res[, i] <- fn_CRfact2OnesZeros(data[,i], in_list_label = in_list_label)
  }
  return(res)
}

fn_CRtabular2indiv <- function(tabular){
  indiv <- tabular[unlist(sapply(1:NROW(tabular), FUN = function(x)rep(x, tabular$Freq[x]))),1:(NCOL(tabular) - 1)]
  return(indiv)
}

fn_CRmatrix2dataframe <- function(data, in_list_label = 'Yes', not_in_list_label = 'No', tabulated = FALSE){
  #
  if(!is.matrix(data)) stop('Object "data" is not a matrix')
  res <- as.data.frame(data)
  lastindx <- NCOL(data)
  if(!tabulated){
    cols <- 1:lastindx
  } else {
    cols <- 1:(lastindx - 1)
    if (!is.numeric(res[,lastindx])) stop('Tabulation column is not numeric')
    res[,lastindx] <- data[,lastindx]
    names(res)[lastindx] <- 'Freq'
  }
  if (!setequal(unique(data[,cols]), c(0,1)) ) stop("Elements of matrix are not exclusively zeros and ones.")
  for (i in cols){
    res[, i] <- factor(data[, i], levels = c(0, 1), labels = c(not_in_list_label, in_list_label) )
  }
  return(res)
}