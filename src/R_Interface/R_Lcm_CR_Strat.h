/* 
 * Copyright (C) 2007-2019 Daniel Manrique-Vallier
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef _R_LCM_CR_STRAT_H_
#define _R_LCM_CR_STRAT_H_

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Utils.h"

extern "C"
{
SEXP R_Create_LCM_CR_Strat(
	SEXP x_flat, SEXP n_strat_plus_NA, SEXP J_plus_strat,
	SEXP n, SEXP K,
	SEXP a_alpha, SEXP b_alpha, SEXP a_lambda, SEXP b_lambda, 
	SEXP min_n, SEXP min_lists, SEXP N_Max_factor,
	SEXP len_buffer, SEXP subsamp
);
}

#endif
