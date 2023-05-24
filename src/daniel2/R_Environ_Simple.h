/* 
 * Copyright (C) 2007-2023 Daniel Manrique-Vallier
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

#ifndef _R_ENVIRON_SIMPLE_H
#define _R_ENVIRON_SIMPLE_H


#include "R.h"
#include "Rinternals.h"
#include "R_ext/Utils.h"
extern "C"{
void finalizer_Env(SEXP ext);
SEXP R_Init_Model(SEXP p);
SEXP R_Update_Model(SEXP p, SEXP int_iter);
SEXP R_Get_Iteration(SEXP p);
SEXP R_Set_Trace(SEXP p, SEXP trace_name);
SEXP R_Get_Param(SEXP p, SEXP param_name);
SEXP R_Set_Param(SEXP p, SEXP param_name, SEXP data);
SEXP R_Get_Trace_List(SEXP p);
SEXP R_Get_Param_List(SEXP p);

SEXP R_Get_Status(SEXP p);
SEXP R_Debug_Trace(SEXP p, SEXP trace);
SEXP R_Get_Trace(SEXP p, SEXP trace);
SEXP R_Reset_Traces(SEXP p);
SEXP R_Activate_Tracing(SEXP p);
SEXP R_Deactivate_Tracing(SEXP p);
SEXP R_Get_Trace_Size(SEXP p);
SEXP R_Change_Trace_Size(SEXP p, SEXP siz);
SEXP R_Change_SubSamp(SEXP p, SEXP subsam);
SEXP R_Set_Seed(SEXP p, SEXP seed);
SEXP R_Activate_Updating_Output(SEXP p);
SEXP R_Deactivate_Updating_Output(SEXP p);
SEXP R_Activate_Chain_Messages(SEXP p);
SEXP R_Deactivate_Chain_Messages(SEXP p);


//MISCELANEOUS FUNCTIONS.
SEXP R_Partial_Contingency_Table(SEXP dataIJ_flat, SEXP levelsJ);
}


//Internal functions (Not to be called from R)
#include "CData.h"
#include "CVariable_Container.h"
#include "CParams_generic.h"
#include "Model_Environ.h"

//Data management functions:
// These functions take arrays in R order (leftmost varying fastest)
//void R_Load_CData(CData& d, SEXP data_list);
//CVariable_Container* R_Array2CVariable_Container(SEXP r_data, const std::string& name);
//void R_Allocate_And_Load_CVariable (SEXP r_data, CVariable_Container& var, const std::string& name);
//void R_Load_CParams_generic(SEXP list, CParams_generic& par);
//void R_Load_CData(CData& d, SEXP data_list);

//More functions
CModel_Environ_Simple_base* get_env(SEXP p);

//maybe a function for creating a generic enviroment object?

//inline SEXP R_Create_Environement_R_Pointer(){
//	SEXP ext = PROTECT(R_MakeExternalPtr((void*)e, R_NilValue, R_NilValue));
//    R_RegisterCFinalizerEx(ext, finalizer_Env, TRUE);
//    UNPROTECT(1);
//	return ext;
//}
// Have to define creators for the environment objects. DO NOT FORGET TO REGISTER THE FINALIZER.
// example:
//SEXP R_Create_Dis_NP_LCM(SEXP x_flat, SEXP J, SEXP n, SEXP levels,SEXP K, SEXP a_alpha, 
//						  SEXP b_alpha, SEXP popsize, SEXP len_buffer, SEXP subsamp){
//	CData_DM* d = new CData_DM();
//	d->Set_Manually(INTEGER(x_flat), *INTEGER(J), *INTEGER(n),INTEGER(levels));
//	env_dis_NP_LCM_t* e = new env_dis_NP_LCM_t(d, new CParams_NP_LCM(d, *INTEGER(K), *REAL(a_alpha), *REAL(b_alpha),  "", ""), 
//		*INTEGER(popsize), *INTEGER(len_buffer), *INTEGER(subsamp), true);
//	SEXP ext = PROTECT(R_MakeExternalPtr((void*)e, R_NilValue, R_NilValue));
//  R_RegisterCFinalizerEx(ext, finalizer_Env, TRUE);
//  UNPROTECT(1);
//	return(ext);
//}
#endif

