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

#include "registration.h"
#include "R_Lcm_CR.h"
#include "R_Lcm_CR_Strat.h"
#include "daniel2/R_Environ_Simple.h"


/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/


static const R_CallMethodDef CallEntries[] = {
    {"R_Activate_Chain_Messages",    (DL_FUNC) &R_Activate_Chain_Messages,    1},
    {"R_Activate_Tracing",           (DL_FUNC) &R_Activate_Tracing,           1},
    {"R_Activate_Updating_Output",   (DL_FUNC) &R_Activate_Updating_Output,   1},
    {"R_Change_SubSamp",             (DL_FUNC) &R_Change_SubSamp,             2},
    {"R_Change_Trace_Size",          (DL_FUNC) &R_Change_Trace_Size,          2},
    {"R_Create_LCM_CR_Basic",        (DL_FUNC) &R_Create_LCM_CR_Basic,        8},
    {"R_Create_LCM_CR_Strat",        (DL_FUNC) &R_Create_LCM_CR_Strat,        14},
    {"R_Deactivate_Chain_Messages",  (DL_FUNC) &R_Deactivate_Chain_Messages,  1},
    {"R_Deactivate_Tracing",         (DL_FUNC) &R_Deactivate_Tracing,         1},
    {"R_Deactivate_Updating_Output", (DL_FUNC) &R_Deactivate_Updating_Output, 1},
    {"R_Get_Iteration",              (DL_FUNC) &R_Get_Iteration,              1},
    {"R_Get_Param",                  (DL_FUNC) &R_Get_Param,                  2},
    {"R_Get_Param_List",             (DL_FUNC) &R_Get_Param_List,             1},
    {"R_Get_Status",                 (DL_FUNC) &R_Get_Status,                 1},
    {"R_Get_Trace",                  (DL_FUNC) &R_Get_Trace,                  2},
    {"R_Get_Trace_List",             (DL_FUNC) &R_Get_Trace_List,             1},
    {"R_Get_Trace_Size",             (DL_FUNC) &R_Get_Trace_Size,             1},
    {"R_Init_Model",                 (DL_FUNC) &R_Init_Model,                 1},
    {"R_Reset_Traces",               (DL_FUNC) &R_Reset_Traces,               1},
    {"R_Set_Seed",                   (DL_FUNC) &R_Set_Seed,                   2},
    {"R_Set_Trace",                  (DL_FUNC) &R_Set_Trace,                  2},
    {"R_Update_Model",               (DL_FUNC) &R_Update_Model,               2},
    {NULL, NULL, 0}
};

extern"C"
void R_init_LCMCR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
