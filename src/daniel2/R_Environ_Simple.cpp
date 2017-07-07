/* 
 * Copyright (C) 2007-2016 Daniel Manrique-Vallier
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

#include "Model_Environ.h"
#include "R_Environ_Simple.h"
#include "dan_math.h"
#include "dan_array_utils.h"
#include <map>
#include <vector>

///////////------------------------------------------------
///////	All these functions will work for any object derived from
/////    CModel_Environ_Simple_base. Reuse at will!
///////////------------------------------------------------
CModel_Environ_Simple_base* get_env(SEXP p) {
	return(reinterpret_cast<CModel_Environ_Simple_base*>(R_ExternalPtrAddr(p)));
}

void finalizer_Env(SEXP ext){
    if (NULL == R_ExternalPtrAddr(ext)) return;
    CModel_Environ_Simple_base* ptr = get_env(ext);//reinterpret_cast<CModel_Environ_Simple_base*>(R_ExternalPtrAddr(ext));
    delete ptr;
    R_ClearExternalPtr(ext);
}


SEXP R_Set_Seed(SEXP p, SEXP seed){
	CModel_Environ_Simple_base* m = get_env(p);
	unsigned int sd = *INTEGER(seed);
	try{
		m->get_CChain().reseed_rng(sd);
	} catch (const std::exception& e){
		DAN_ERR_NOABORT("Cannot set rng seed to %d", sd);
	}
	return(p);
}
SEXP R_Activate_Updating_Output(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->activate_updating_output();
	return(p);
}
SEXP R_Deactivate_Updating_Output(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->deactivate_updating_output();
	return(p);
}
SEXP R_Activate_Chain_Messages(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->get_CChain().messages_on();
	return(p);
}
SEXP R_Deactivate_Chain_Messages(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->get_CChain().messages_off();
	return(p);
}


SEXP R_Init_Model(SEXP p){
	CModel_Environ_Simple_base* m;
	try {
		m = get_env(p);
		m->init_model();
	} catch (const std::exception& e) {
		std::string s = "Model initialization failed (";
		s += e.what();
		s += ")";
		DAN_ERR_NOABORT(s.c_str());
		p = R_NilValue;
	} 
	return(p);
}
SEXP R_Update_Model(SEXP p, SEXP int_iter){
	CModel_Environ_Simple_base* m = get_env(p);
	try{
		m->update_model(*INTEGER(int_iter));
	} catch (const std::exception& e){
		DAN_ERR_EXIT("%s", e.what());
	}
	return(p);
}

SEXP R_Get_Iteration(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	SEXP iter;
	PROTECT(iter = allocVector(INTSXP,1));
	*INTEGER(iter) = m->get_iteration();
	UNPROTECT(1);
	return(iter);
}

SEXP R_Get_Status(SEXP p){
	//<<<---- TODO: FINISH THIS FUNCTION.
	//return a vector with the statust of the chain
	//Returns a vector c(<iteration>, <initialized>, <buffer size>, <buffer used>, <tracer activated>, <thinning>) 
	CModel_Environ_Simple_base* m = get_env(p);
	SEXP status;
	PROTECT(status = allocVector(INTSXP,6));
	INTEGER(status)[0] = m->get_iteration();
	INTEGER(status)[1] = (m->get_model_state() == CModel_Environ_Simple_base::MODEL_INITIALIZED ||  m->get_model_state() == CModel_Environ_Simple_base::UPDATING) ? 1 : 0;
	INTEGER(status)[2] = m->Tracer().get_capacity();
	INTEGER(status)[3] = m->Tracer().get_size();
	INTEGER(status)[4] = m->get_tracing() ? 1 : 0;
	INTEGER(status)[5] = m->get_thinning();
	UNPROTECT(1);
	return(status);
}


SEXP R_Set_Trace(SEXP p, SEXP trace_name){
	CModel_Environ_Simple_base* m = get_env(p);
	char* name;
	name = const_cast<char*>(CHAR(STRING_ELT(trace_name,0)));
	try{
		m->Tracer().set_trace(name);
	} catch (const std::exception& e){
		DAN_ERR_NOABORT(e.what());
	}
	return(p);
}

SEXP R_Debug_Trace(SEXP p, SEXP trace){
	CModel_Environ_Simple_base* m = get_env(p);
	char* key = const_cast<char*>(CHAR(STRING_ELT(trace,0)));

	std::vector<int> dims = m->Tracer().get_dims(key);
	int ndims = dims.size();
	int nelems = m->Tracer().get_size_elems(key);

	DAN_PRINTF("Size elements = %d\n", nelems);
	DAN_PRINTF("Size bytes = %d\n", m->Tracer().get_size_bytes(key));
	DAN_PRINTF("ndims = %d\n", ndims);
	DAN_PRINTF("dims : ");
	for (int i = 0; i < ndims; i++){
		DAN_PRINTF("%d ", dims[i]);
	}
	DAN_PRINTF("\n");
	//Test one of the traces

	return(p);
}


SEXP R_Get_Trace(SEXP p, SEXP trace){
	CModel_Environ_Simple_base* m = get_env(p);
	SEXP r = R_NilValue; //int or double
	SEXP ret_int_dims = R_NilValue; //int
	char* key = const_cast<char*>(CHAR(STRING_ELT(trace,0)));
	if (!m->Tracer().check_trace(key))
		return(R_NilValue);
	std::vector<int> dims = m->Tracer().get_dims(key); //NOTE THAT FIRST DIMENSION IS THE *CAPACITY*, NOT ACTUAL SIZE
	int ndims = dims.size();
	int nelems = m->Tracer().get_size_elems(key);
	CPar_Data_Type& type = m->Tracer().get_data_type(key);
	switch (type.get_data_type()){
		case CPar_Data_Type::T_DOUBLE: 
			PROTECT(r = allocVector(REALSXP, nelems));
			m->Tracer().Copy_trace(key, REAL(r));
			break;
		case CPar_Data_Type::T_INT:
			PROTECT(r = allocVector(INTSXP, nelems));
			m->Tracer().Copy_trace(key, INTEGER(r));
			break;
		default:
			//other types not implemented at the moment.
			DAN_ERR_EXIT("Data type not implemented\n");
			break;
	}
	// Note that the vector is in rightmost-varying-fastest. We set the "dim"
	// vector in reverse order and use "aperm(x, ndims:1)"
	PROTECT(ret_int_dims = allocVector(INTSXP, ndims));
	for (int i = 1; i < ndims; i++){
		INTEGER(ret_int_dims)[ndims - i -1] = dims[i];
	}
	INTEGER(ret_int_dims)[ndims - 1] = m->Tracer().get_size();//->get_traces_current_size();
	setAttrib(r, install("dim"), ret_int_dims);
	UNPROTECT(2);
	return(r);
}

SEXP R_Get_Param(SEXP p, SEXP key_){
	CModel_Environ_Simple_base* m = get_env(p);
	SEXP r = R_NilValue; //int or double
	SEXP ret_int_dims = R_NilValue; //int
	char* key = const_cast<char*>(CHAR(STRING_ELT(key_,0)));
	if(!m->check_param_key(key)) return ( R_NilValue );
	CPar_defs& par = m->get_param_container(key); 
	const std::vector<int>& dims = par.get_dim_lengths();
	int ndims = dims.size();
	switch (par.get_data_type().get_data_type()){ 
		case CPar_Data_Type::T_DOUBLE:
			PROTECT(r = allocVector(REALSXP, par.get_size_elems()));
			par.copy_raw_data(REAL(r));
			break;
		case CPar_Data_Type::T_INT:
			PROTECT(r = allocVector(INTSXP, par.get_size_elems()));
			par.copy_raw_data(INTEGER(r));
			break;
		default:
			DAN_ERR_NOABORT("Can't get variable. Not implemented data type. \n");
			return(R_NilValue);
	}
	PROTECT(ret_int_dims = allocVector(INTSXP, ndims));
	for (int i = 0; i < ndims; i++){
		INTEGER(ret_int_dims)[ndims - i - 1] = dims[i];
	}
	setAttrib(r, install("dim"), ret_int_dims);
	UNPROTECT(2);
	return(r);
}

 
//<<--- TEST THIS FUNCTION.
SEXP R_Set_Param(SEXP p, SEXP param_name, SEXP r_data){
	//Copy provided data r_data to parameter "param_name". Verifies
	// data type and dimensions.
	CModel_Environ_Simple_base* m = get_env(p);
	//SEXP ret_int_dims = R_NilValue; //int. Dimensions of the array
	//Get and verify parameter container.
	char* key = const_cast<char*>(CHAR(STRING_ELT(param_name,0)));
	if(!m->check_param_key(key)){
		DAN_ERR_EXIT("No such parameter (%s)\n", key);
	}
	CPar_defs& par = m->get_param_container(key);
	//get and verify input data type
	CPar_Data_Type::data_type_t type;
	void* raw_pointer;
	std::string err;
	switch( TYPEOF(r_data) ){
	case INTSXP:
		type = CPar_Data_Type::T_INT;
		raw_pointer = (void*) INTEGER(r_data);
		break;
	case REALSXP:
		type = CPar_Data_Type::T_DOUBLE;
		raw_pointer = (void *) REAL(r_data);
		break;
	default:
		err = "Not implemented data type reader for " + par.get_name();
		//throw std::runtime_error(err);
		DAN_ERR_EXIT("%s\n", err.c_str());
		break;
	}
	if (type != par.get_data_type().get_data_type()){
		DAN_ERR_EXIT("Incorrect data type for parameter %s (should be %s)\n", 
			par.get_name().c_str(),
			par.get_data_type().get_type_name().c_str()
		);
	}
	//Get and verify dimensions
	int size_elems = Rf_length(r_data);
	SEXP r_dims = getAttrib(r_data, R_DimSymbol);
	int n_r_dims = Rf_length(r_dims);
	n_r_dims = n_r_dims == 0 ? 1 : n_r_dims; //in case the input is just a vector ("0" dims)
	int n_c_dims = par.get_dims();
	std::vector<int> dims(n_r_dims);
	std::copy_backward(INTEGER(r_dims), INTEGER(r_dims) + n_r_dims, dims.begin());
	if (n_r_dims != n_c_dims){
		DAN_ERR_EXIT("Incorrect number of dimensions (provided %d; should be %d)\n", n_r_dims, n_c_dims);
	}
	if (n_r_dims > 1 && !std::equal(dims.begin(), dims.end(), par.get_dim_lengths().begin())){
		DAN_ERR_EXIT("Dimensions don't match.\n");
	}
	if (size_elems != par.get_size_elems()){
		//Last verification: the number of elements (and bytes!) has to be the same.
		DAN_ERR_EXIT("Incorrect number of elements (provided %d; should be %d).\n", size_elems, par.get_size_elems());
	}
	//copy data. Transposition between R and C here.
	dan_transpose_untyped(raw_pointer, par.get_data_base(), 
		par.get_dims(), &dims.front(), 
		par.get_size_dataelem()
	);
	return(R_NilValue);
}


SEXP R_Get_Trace_List(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	const std::vector<std::string> &t =  m->Tracer().get_trace_keys();//get_trace_keys();
	int n = t.size();
	SEXP names;
	PROTECT(names = allocVector(STRSXP,n));
	for(int i = 0; i < n; i++)   
		SET_STRING_ELT(names, i, mkChar(t[i].c_str())); 
	UNPROTECT(1);
	return(names);
}

SEXP R_Get_Param_List(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	const std::vector<std::string> &t =  m->get_param_keys();
	int n = t.size();
	SEXP names;
	PROTECT(names = allocVector(STRSXP,n));
	for(int i = 0; i < n; i++)   
		SET_STRING_ELT(names, i, mkChar(t[i].c_str())); 
	UNPROTECT(1);
	return(names);
}

SEXP R_Reset_Traces(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->Tracer().Clear();//reset_traces();
	return(p);
}

SEXP R_Activate_Tracing(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->activate_tracing();
	return(p);
}

SEXP R_Deactivate_Tracing(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	m->deactivate_tracing();
	return(p);
}

SEXP R_Get_Trace_Size(SEXP p){
	CModel_Environ_Simple_base* m = get_env(p);
	SEXP r;
	PROTECT(r = allocVector(INTSXP, 1));
	*INTEGER(r)= m->Tracer().get_size();//get_traces_current_size();
	UNPROTECT(1);
	return(r);
}

SEXP R_Change_Trace_Size(SEXP p, SEXP siz){
	CModel_Environ_Simple_base* m = get_env(p);
	try {
		m->recreate_trace(*INTEGER(siz));
	} catch (const std::exception& e){
		DAN_ERR_NOABORT(e.what());
	}
	return(R_NilValue);
}

SEXP R_Change_SubSamp(SEXP p, SEXP subsam){
	CModel_Environ_Simple_base* m = get_env(p);
	m->reset_subsamp(*INTEGER(subsam));
	return(R_NilValue);
}

SEXP R_Partial_Contingency_Table(SEXP dataIJ_flat, SEXP levelsJ){
	int J = length(levelsJ);
	int JN = length(dataIJ_flat);
	//int N = JN / J;
	int* data_flat = INTEGER(dataIJ_flat);
	typedef std::map<std::vector<int>, int> tabla;
	tabla s;
	for (int* it = data_flat; it != data_flat + JN; it += J){
		++s[std::vector<int>(it, it + J)];
	}
	//prepare the return list
	int M = s.size();
	SEXP ans, cells, Freq, dims;
	PROTECT(cells = allocVector(INTSXP, J * M));
	PROTECT(Freq = allocVector(INTSXP, M));
	PROTECT(dims = allocVector(INTSXP, 2));
	PROTECT(ans = allocVector(VECSXP, 2));
	int m = 0;
	int* cell = INTEGER(cells);
	for (tabla::iterator it = s.begin(); it != s.end(); ++it, ++m, cell += J){
		std::copy(&(it->first[0]), &(it->first[J]), cell);
		INTEGER(Freq)[m] = it->second;
	}
	INTEGER(dims)[0] = J; //COL major order.
	INTEGER(dims)[1] = M;
	setAttrib(cells, install("dim"), dims);
	SET_VECTOR_ELT(ans, 0, cells);
	SET_VECTOR_ELT(ans, 1, Freq);
	UNPROTECT(4);
	return(ans);
}


CVariable_Container* R_Array2CVariable_Container(SEXP r_data, const std::string& name){
	//INTSXP and REALSXPs into a new CVariable_Container
	//If exception, the object is not created.
	CPar_Data_Type::data_type_t type;
	void* raw_pointer;
	std::string err;
	switch( TYPEOF(r_data) ){
	case INTSXP:
		type = CPar_Data_Type::T_INT;
		raw_pointer = (void*) INTEGER(r_data);
		break;
	case REALSXP:
		type = CPar_Data_Type::T_DOUBLE;
		raw_pointer = (void *) REAL(r_data);
		break;
	default:
		err = "Not implemented data type reader for " + name;
		throw std::runtime_error(err);
		break;
	}
	//Determine dimensions.
	int size_elems = Rf_length(r_data);
	SEXP r_dims = getAttrib(r_data, R_DimSymbol);
	int n_dims = Rf_length(r_dims);
	std::vector<int> v(n_dims);
	if(n_dims == 0) {
		n_dims = 1;
		v.push_back(size_elems);
	} else {
		int* lengths = INTEGER(r_dims);
		std::copy(lengths, lengths + n_dims, v.begin());
	}
	//allocate and copy data. -> exception thrown if failed.
	CVariable_Container* var = new CVariable_Container(type, n_dims, name);

	try{
		var->allocate_space(v);
	} catch (const std::exception& e){
		err = "Cannot allocate space for " + name + ". Object destroyed.";
		delete var;
		throw std::runtime_error(err);
	}
	std::reverse(v.begin(), v.end());
	dan_transpose_untyped(raw_pointer, var->get_data_base(), 
		var->get_dims(), &v.front(), 
		var->get_size_dataelem()
	);
	return var;
}

void R_Allocate_And_Load_CVariable (SEXP r_data, CVariable_Container& var, const std::string& name){
	//reads INTSXP and REALSXPs into CVariable_Container
	void* raw_pointer;
	std::string err;
	//determine data type
	switch( TYPEOF(r_data) ){
	case INTSXP:
		if (var.get_data_type().get_data_type() != CPar_Data_Type::T_INT){
			//abort. Bad data type
			err = name + ": wrong data type (tried to load INT)";
			throw std::runtime_error(err);
		}
		raw_pointer = (void*) INTEGER(r_data);
		break;
	case REALSXP:
		if (var.get_data_type().get_data_type() != CPar_Data_Type::T_DOUBLE){
			err = name + ": wrong data type (tried to load DOUBLE)";
			throw std::runtime_error(err);
		}
		raw_pointer = (void *) REAL(r_data);
		break;
	default:
		throw std::runtime_error("Not implemented data type reader");
		break;
	}
	//Determine dimensions.
	int size_elems = Rf_length(r_data);
	SEXP r_dims = getAttrib(r_data, R_DimSymbol);
	int n_dims = Rf_length(r_dims);
	std::vector<int> v(n_dims);
	if(n_dims == 0) {
		n_dims = 1;
		v.push_back(size_elems);
	} else {
		int* lengths = INTEGER(r_dims);
		std::copy(lengths, lengths + n_dims, v.begin());
	}

	//allocate and copy data.
	//WE TRANSPOSE IT WHILE COPYING!
	var.allocate_space(v);
	std::reverse(v.begin(), v.end());
	dan_transpose_untyped(raw_pointer, var.get_data_base(), 
		var.get_dims(), &v.front(), 
		var.get_size_dataelem()
	);
}

void R_Load_CParams_generic(SEXP list, CParams_generic& par){
	if ( TYPEOF(list) != VECSXP ){
		throw std::runtime_error("Bad R data type (should be list)");
	}
	int n_variables = Rf_length(list); //get list length
	SEXP names = getAttrib(list, R_NamesSymbol); 
	for (int i = 0; i < n_variables; i++){
		SEXP r_name = STRING_ELT(names, i);
		const char* name = CHAR(r_name);
		if (Rf_StringBlank(r_name)){
			throw std::runtime_error("Not named element in list");
		}
		SEXP data_object = VECTOR_ELT(list, i);
		//throws an exception if error.
		CVariable_Container* var = R_Array2CVariable_Container(data_object, name);
		par.add(*var);
	}
}

void R_Load_CData(CData& d, SEXP data_list){
	int n_variables = Rf_length(data_list); //get list length
	// - Obtain list of names (R - vector of strings: STRING_ELT)
	SEXP names = getAttrib(data_list, R_NamesSymbol); 
	//Read a list
	for (int i = 0; i < n_variables; i++){
		const char* name = CHAR(STRING_ELT(names, i));
		SEXP data_object = VECTOR_ELT(data_list, i);
		if ( !d.get_data_container().check_key(name) ) break;
		//Load the data
		R_Allocate_And_Load_CVariable(data_object, d[name], name);
		d._Mark_Var_As_Read(name); //<-important. Keeps the status of the CData object consistent.
	}
	if ( !d.is_loaded() ){
		throw std::runtime_error("Unfinished loading in CData object. Object might be inconsistent.");
	}
}

