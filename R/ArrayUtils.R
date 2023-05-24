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

###############################################################################
#Functions for manipulating matrices and dataframes with categorical values
# (c) Daniel Manrique-Vallier 2012
###############################################################################

fn_df_nlevels <- function(d){
    sapply(names(d), FUN = function(x)nlevels(d[,x]))
}

fn_df_discretize <- function(d, cols = 1:NCOL(d)){
  for (i in cols){
    d[,i] <- factor(d[,i])
  }
  return(d)
}
fn_df_levels <- function(d, cols = 1:NCOL(d)){

  res <- list()
  for (i in cols){
    res[[names(d)[i]]] <- levels(d[1,i])
  }
  return(res)
}

fn_numeric_to_factor <- function(
  x, list_levls = rep(list(1:2), NCOL(x)), cols = 1:NCOL(x), missing = -1, offset = 0){
  x <- as.data.frame(x)
  for (c in cols){
    lv <- list_levls[[c]]
    nlv <- length(lv)
    x[,c] <- factor(x[,c], levels = 1:nlv + offset)
    levels(x[,c]) <- lv
  }
  return(x)
}

fn_apply_levels_from <- function(dest, from, cols = 1:NCOL(from)){
  dest <- as.data.frame(dest)
  lev <- sapply(cols, function(x)levels(from[1,x]), simplify=FALSE)
  for (i in cols){
    dest[,i] <- factor(dest[,i], levels = lev[[i]])
  }
  names(dest[,cols]) <- names(from[,cols])
  return(dest)
}

fn_dataframe2num_matrix <- function(d, offset = -1, missing = -1, C_style= FALSE){
  offset <- as.integer(offset); missing <- as.integer(missing)
  r <- data.matrix(d) + offset
  r[is.na(r)] <- missing
  if(C_style){
    r <- t(r)
  }
  class(r) <- c('matrix', 'DMVmatrix')
  attr(r, 'levels') <- fn_df_levels(d)
  attr(r, 'na.code') <- missing
  attr(r, 'offset') <- offset
  attr(r, 'C_style') <- C_style
  return(r)
}

fn_df_2_C_Data <- function(df, prefix = 'array', filename = 'console'){
  #Use this function for generating test data. Works with arrays if factors.
  breakdown <- fn_dataframe2num_matrix(df, offset = -1, missing=-1, C_style=T)
  n <- NROW(df)
  J <- NCOL(df)
  l1 <- paste('int ',prefix,'_data_raw[]= {', paste(breakdown, collapse=', '), '};\n', sep='')
  l2 <- paste('int ', prefix,'_levels[]= {', paste(fn_df_nlevels(df), collapse=', '), '};\n', sep ='')
  l3 <- paste('int ',prefix, '_n_glob = ', n,';\n', 'int ', prefix,'_J_glob=', J,';\n', sep='')
  if(filename != 'console'){
    fl = file(description = filename, open = 'w')
    writeLines(strwrap(paste(l1, l2, l3, sep ='\n')), con = fl)
    close(fl)
  } else {
    writeLines(strwrap(paste(l1, l2, l3, sep ='\n')))
  }
}
