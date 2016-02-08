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

#ifndef _CVARIABLE_CONTAINER_H_
#define _CVARIABLE_CONTAINER_H_

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "CPar_Data_Type.h"
//#include "daniel2/dan_sys.h"

class CVariable_Container{
public:
	//STATES OF THE OBJECT:
	// UNINITIALIZED -> Only name set
	// METADATA_LOADED -> Name + type + dims loaded.
	// ALLOCATED -> Everything set. Object ready to use.
	typedef enum {UNINITIALIZED, METADATA_LOADED, ALLOCATED} state_t; 
	state_t get_state(){ return state; }
	//Data type
	CPar_Data_Type data_type;

	// constructors:
	//1) with external storage
	CVariable_Container(CPar_Data_Type::data_type_t type, void* existing, const std::string& sname) : data_type(type), state(UNINITIALIZED) {
		construct_class(sname);
		add_existing_scalar(existing);
	}
	//2) with internal storage.
	CVariable_Container(CPar_Data_Type::data_type_t type, const std::vector<int>& lengths, const std::string& sname) 
		:  data_type(type), state(UNINITIALIZED) {
		construct_class(sname);
		alloc(lengths);
	}
	//3) without allocation. Just name it, type it, and load the number of dimensions.
	CVariable_Container(CPar_Data_Type::data_type_t type, int dims, const std::string& sname) 
		:  data_type(type), state(UNINITIALIZED) {
		construct_class(sname);
		this->dims = dims;
		state = METADATA_LOADED;
	}
	//4) Copy constructor. Note that this constructor ignores the "existing" flag. It just allocates storage.
		// we can change this later using the "swap_internal2existing_scalar".
	CVariable_Container(CVariable_Container& orig) : data_type(orig.get_data_type().get_data_type()), state(UNINITIALIZED){
		construct_class(orig.name);
		if (orig.state != ALLOCATED) return; //If the state is not allocated in the original, don't do it in the copy.
		this->alloc(orig.dim_lengths);
		this->copy_from_raw(orig.data_base);
	}
	//Destructor. Deallocates space.
	virtual ~CVariable_Container(){
		if (state == ALLOCATED && !this->existing) {
			operator delete(data_base);
			if (this->dims >1) operator delete(data);
		}
	}

	//accessors. Read only.
	//1) Properties that are always valid
	CPar_Data_Type& get_data_type(){ return data_type; }
	int	get_dims()		{ return dims; }
	int get_size_dataelem(){ return data_type.get_bytes(); }
	std::string get_name(){	return name; }
	int	get_bytes_per_elem() { return data_type.get_bytes();}
	
	//2) Properties valid only after allocation. 
	std::vector<int>& get_dim_lengths();
	void*	get_data();	
	void*	get_data_base();
	int		get_size_bytes();
	int		get_size_elems();
	bool	get_existing();	

	//Actions
	void swap_internal2external_scalar(void* existing_scalar);
	void copy_from_raw(void* src);
	void copy_raw_data(void* dst);
	template <typename _T>
	void fill(const _T& value){
		//minimal type check: verifies that the object's datatype length is at least a multiple of _T.
		if (data_type.get_bytes() % sizeof(_T) != 0) throw std::runtime_error("Types are not compatible");
		int n_copies = size_bytes / sizeof(_T);
		std::fill( (_T*)data_base, (_T*)data_base + n_copies, value);
	}
	//For using in objects int the "METADATA_LOADED" state
	void allocate_space(const std::vector<int>& lengths);
	void allocate_space(int dims, ...);
	void register_data(const std::vector<int>& lengths, void* raw_data, bool cleanup = true);

	//A multimensional-to-flat address converter. Can give more efficient memory access when not in a linear sequence.
	template<typename _T>
	_T& arr_elem(int& i){
		return ( (reinterpret_cast<_T*>(data_base))[i] );
	}
	template<typename _T>
	_T& arr_elem(int& i, int& j){
		return ( (reinterpret_cast<_T*>(data_base))[i * dim_lengths[1] + j] );
	}
	template<typename _T>
	_T& arr_elem(int& i, int& j, int& k){
		return ( (reinterpret_cast<_T*>(data_base))[(i * dim_lengths[1] + j)*dim_lengths[2] + k] );
	}
	template<typename _T>
	_T& arr_elem(int& i, int& j, int& k, int& l){
		return ( (reinterpret_cast<_T*>(data_base))[((i * dim_lengths[1] + j)*dim_lengths[2] + k)*dim_lengths[3] + l] );
	}
	template<typename _T>
	_T& arr_elem(int& i, int& j, int& k, int& l, int& m){
		return ( (reinterpret_cast<_T*>(data_base))[(((i * dim_lengths[1] + j)*dim_lengths[2] + k)*dim_lengths[3] + l)*dim_lengths[4] + m] );
	}
	////debugging methods
	//void dump_content(){
	//	for (int i = 0; i < this->size_elems; i++){
	//		switch(data_type.get_data_type()){
	//		case CPar_Data_Type::T_DOUBLE:
	//			DAN_PRINTF("%.3f ", ((double*)this->data_base)[i]);
	//			break;
	//		case CPar_Data_Type::T_INT:
	//			DAN_PRINTF("%d ", ((int*)this->data_base)[i]);
	//			break;
	//		}
	//	}
	//}
private:
	//properties
	state_t			state; //state of the object.
	std::string		name;
	int				dims; //number of dimensions
	std::vector<int>dim_lengths; //size of each dimension.
	void*			data; //pointer to the data structure (THE ARRAY OF POINTERS). HAS TO BE RECASTED IN THE CALLING CLASS.
	void*			data_base; //the base address of the data (the actual memory obtained through malloc)
	int				size_bytes; //allocated size (of the actual data) in bytes
	int				size_elems; //number of allocated elements. (size_bytes/sizeof(numeric_type))
	bool			existing; //if the storage was assigned elsewhere.


	/////////////
	//Private Methods
	CVariable_Container();// Disable default constructor. 
	void add_existing_scalar(void* var); //only called from constructor.
	void alloc(const std::vector<int>& lengths);
	void construct_class(const std::string& sname){		
		dims = size_elems  = size_bytes = 0;
		data_base = data = NULL;
		existing = false;
		state = UNINITIALIZED;
		this->name = sname;
	}
	void abort_if_not_allocated();
};

#endif

