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

/*
	CData_CR class declaration.
	This function holds the data for the computations.
	It includes accessors and a function for reading the data from files
*/

#ifndef _CDATA_H
#define _CDATA_H

#include <string>
#include <cstdarg>
#include <stdexcept>
#include "CParams_generic.h"
#include "CVariable_Container.h"

class CData {
	//1) Implement a *PRIVATE* method "Declare" that uses the method _Declare_Variable to
	//		add the variables to the container, but doesn't allocate them. Important: this function is to be called from the constructor,
	//		so keep it private so that derived classes use their own implementation!
	//1.1) The method "Declare" should only initialize the variables corresponding to the class. If the class is derived from other CData child, 
	//		the parent's declaration will be called from its own constructor.
	//2) Call this method from the *CHILD* class' constructor. 
	//3) Implement a (virtual) function Close_Loading() that calls (PARENT)::Close_Loading() and register all the pointers 
	//		that expose data in the child CData object, computes all the necessary statistics, and load
	//		any other needed variable (e.g. dimensions of the arrays). To get the pointers, use "_get_pointer".
	//4) _Mark_Var_As_Read(key) automatically calls Close_Loading() when the last variable has loaded.
	//5) Create a method called "SetManually( ...pars with data... )" that uses _Load_Variable(...) 
	//		to allocate and fill the variables. 
	//6) At any moment the vector "vars_to_read" (readable through _need_to_load()) indicates which variables
	//		have to be loaded before closing the object for reading. _Load_Variable(..) keeps this updated
	//		automatically. If you use other method, this has to be updated manually with a call to _Mark_Var_As_Read(key).


public:
	CData(){}; //Default constructor. Nothing to do here.
	
	//<---TODO: [ ] 2016-02-09: Check the effect of introducing this destructor.
	virtual ~CData(){}; //No destructor needed because all allocations are placed into the member object "data_container"
	CVariable_Container& operator[] (const std::string& key){
		if ( !data_container.check_key(key) ){
			throw std::runtime_error(key + " is not a registered variable in CData object");
		}
		return data_container[key];
	}
	bool is_loaded() { return vars_to_read.size() == 0; }
	CParams_generic& get_data_container(){ return data_container;}
	void* _Load_Variable(const std::string key, void* raw_data, const std::vector<int> lengths){
		if ( !data_container.check_key(key) ){
			throw std::runtime_error("Variable " + key + " not defined in CData object");
		}
		void* dst_pointer = data_container.alloc_delayed(key, lengths);
		data_container[key].copy_from_raw(raw_data);
		//If the key is in the list of vars to read, we mark it as read.
		if ( std::find(vars_to_read.begin(), vars_to_read.end(), key) != vars_to_read.end()){
			_Mark_Var_As_Read(key);
		}
		return dst_pointer;
	}

	void* _Load_Variable(const std::string& key, void* raw_data, ...){
		int i;
		int dims = data_container[key].get_dims();//.get_dims(key);
		std::vector<int> lengths( dims );
		std::va_list args;
		va_start(args, raw_data);
		for (i = 0; i < dims; i++){
			lengths[i] = va_arg(args, int);
		}
		va_end(args);
		return _Load_Variable(key, raw_data, lengths);
	}

	std::vector<std::string> _need_to_load(){
		return vars_to_read;
	}

	void _Mark_Var_As_Read(const std::string& key){
		if (vars_to_read.size() == 0) {
			throw std::runtime_error("CData already closed.");
		}
		std::vector<std::string>::iterator it = std::find(vars_to_read.begin(), vars_to_read.end(), key);
		if (it != vars_to_read.end()){
			vars_to_read.erase(it);
		}
		if ( vars_to_read.size() == 0 ){
			Close_Loading();
		}
	}

	void* _Substitute_Variable(CVariable_Container& new_var){
		//Be careful with this function.
		//it doesn't automatically update the pointer that the child class might have set to the 
		// the data areas. THOSE HAVE TO BE RE SET MANUALLY.
		const std::string& key = new_var.get_name();
		if (!data_container.check_key(key)){
			throw std::runtime_error("Cannot substitute. Variable doesn't exist");
		}
		CPar_defs& orig = data_container.get_dataObject(key);
		if (orig.get_state() != CPar_defs::ALLOCATED){
			throw std::runtime_error("Cannot substitute. Variable still not allocated");
		}
		data_container.erase_variable(key);
		return(data_container.add(new_var));

	}
protected:
	void* _get_pointer(const std::string& key){
		return data_container[key].get_data();
	}
	void _set_internal_pointer(const std::string& key, void* dst_pointer_address){
		*((void**)dst_pointer_address) = data_container[key].get_data();
	}
	void _Declare_Variable(const std::string& name, CPar_Data_Type::data_type_t type, int dims, bool derived){
		data_container.add_no_alloc(name, type, dims);
		if (!derived){
			vars_to_read.push_back(name);
		}
	}
	void* _Declare_and_Allocate_derived(const std::string& key, CPar_Data_Type::data_type_t type, int dims, ...){
		int i;
		std::vector<int> lengths( dims );
		std::va_list args;
		va_start(args, dims);
		for (i = 0; i < dims; i++){
			lengths[i] = va_arg(args, int);
		}
		va_end(args);
		_Declare_Variable(key, type, dims, true);
		void* dst_pointer = data_container.alloc_delayed(key, lengths);
		return dst_pointer;
	}
	std::vector< std::string > vars_to_read; //variables that need to be read from somewhere (files, memory, who knows...)
	CParams_generic data_container;
	//pure virtual methods (HAVE TO IMPLEMENT IN CHILD CLASS)
	//virtual void Declare() = 0; //declare data structure, but do not allocate. Use "Declare_Variable" to do this.
	virtual void Close_Loading(){
	
	};
};

#endif  //_CDATA_H

///////////////
// EXAMPLE
//////////////
//class CData_DM : public CData{
//public:
//	CData_DM(){
//		Declare();
//	}
//	void Set_Manually(int *x_flat, int J, int n, int *levels){
//		_Load_Variable("x", x_flat, n, J);
//		_Load_Variable("levelsJ", levels, J);
//	}
//	//data storage. These should be private, but this is going to be more efficient.
//	int ** x; //responses vector
//	int *levelsJ; // (levelsJ[j] = #levels of variable j)
//	int L; //max #of levels (for dimentions)
//	int J;
//	int n; //sample size
//protected:
//	void Close_Loading(){
//		CData::Close_Loading();
//		_set_internal_pointer("x", &x);
//		_set_internal_pointer("levelsJ", &levelsJ);
//		J = data_container["x"].get_dim_lengths()[1];
//		n = data_container["x"].get_dim_lengths()[0];
//		L = *std::max_element(levelsJ, levelsJ + J);
//	}
//private:
//	void Declare(){
//		_Declare_Variable("x", CPar_Data_Type::T_INT, 2, false);
//		_Declare_Variable("levelsJ", CPar_Data_Type::T_INT, 1, false);
//	}
//};
