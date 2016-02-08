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

#ifndef _CPARAM_H
#define _CPARAM_H

#include "CParams_generic.h"
#include <string>
#include <gsl/gsl_rng.h>
#include <cstdarg>

// Guidelines:
// - all allocations must be done in the "storage" container. This takes care of cleaning, copying and a lot of cool stuff
// - all local values that define the status of the object must be registered in "storage" using "add_existing_scalar".
// - All the variables that define the state of the object MUST BE REGISTERED in the "storage" container using the "_add"-functions
// - ... using the storage object in this form allows to use the serialization, copying and copy constructing functions without extra implementation.
// - The operator "=" just copies the container. Shouldn't be overriden./
// - Each subclass should implement a (non virtual) private "class_construct" method for handling the allocations and initializations that will be called
//   from the constructors.
// - The destructor is virtual, but one shouldn't need to override it. If all the guidelines were followed, no allocations are done outside "storage".
// - In "class_construct", all allocations into the container have to be included into a "try" block that catches a bad_alloc exception.
//		If an allocation fails, we need to deallocate everything and throw another bad_alloc exception with asuitable message. 
// - If the allocation (in the constructor) of any of the variables fails, all the already allocated variaables
//		are properly destroyed when the member object storage is automatically destroyed.


class CParam {
public:
	//parameters
	double loglike;

	//functions
	CParam() {
		//call class_construct() in child class.
	}; 
	CParam& operator = (CParam  &orig){
		if(&orig == this) return (*this);
		//this assumes that the objects are compatible. This is the programmer's responsibility.
		this->storage = orig.storage;
		return(*this);
	}
	virtual ~CParam(){}; //Destructor
	virtual void initizalize(gsl_rng *r) = 0; //(Initialize the parameters). Implement in child class. 
	CParams_generic& get_container(){return(storage);}
	int get_blob_size(){ return (storage.get_blob_size());}
	void get_blob(unsigned char* buffer) { storage.serialize_out(buffer); }
	void set_blob(unsigned char* buffer){ storage.serialize_in(buffer); }
	CVariable_Container& operator[] (const std::string& key){
		return storage[key];
	}
	void set_parameter(const std::string& key, void *raw_data){
		//storage.fill_from_raw(key, raw_data);
		storage[key].copy_from_raw(raw_data);
	}
	std::vector<std::string> get_var_list(){
		return storage.get_registered_keys();
	}
	bool check_key(const std::string& key){
		return storage.check_key(key);
	}
protected:
	CParams_generic storage;
	void _add_parameter(const std::string& key, const CPar_Data_Type::data_type_t& type, void* dst_variable_address, int dims, ...){
				int i;
		std::vector<int> lengths( dims );
		va_list args;
		va_start(args, dims);
		for (i = 0; i < dims; i++){
			lengths[i] = va_arg(args, int);
		}
		va_end(args);
		*((void**)dst_variable_address) = storage.add(key, type, lengths);
	}
	void _add_existing_scalar(const std::string& key, const CPar_Data_Type::data_type_t& type, void* existing_address){
		storage.add_existing_scalar(key, type, existing_address);	
	}
	void _remove_parameter(const std::string& key){
		//Be careful with pointers that point to these data!!!
		storage.erase_variable(key);
	}
private:
	//Create a class_construct() method and call it from the child constructor.
};


#endif
