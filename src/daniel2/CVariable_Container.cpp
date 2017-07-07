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

#include <cstdlib>
#include <cstdarg> 
#include "dan_array_utils.h"
#include "CParams_generic.h"
#include <iostream>
#include <iterator>
#include <stdexcept>



void CVariable_Container::alloc(const std::vector<int>& lengths){
	//allocates space and adds new key.
	if (lengths.size() < 1 || state == ALLOCATED) return;
	int size_elem = data_type.get_bytes();
	this->size_bytes = data_type.get_bytes(); 
	for (unsigned int i = 0; i < lengths.size(); i++) size_bytes *= lengths[i];
	this->data_base = ::operator new (size_bytes); //This throws bad_alloc if fails.

	//Exception handling: if the first allocation works, but the second fails, we need to deallocate
	// the first.
	try{
		this->data = dan_flat2arrayND_cpp(data_base, size_elem, lengths);
	} catch (const std::bad_alloc& e){
		//cleanup: delete main allocation
		::operator delete( data_base );
		throw std::bad_alloc(e);
	}

	//this->size_of_elem = size_elem;
	this->dims = lengths.size();
	this->existing = false;
	this->size_elems = size_bytes / size_elem;
	std::copy(lengths.begin(), lengths.end(), std::back_inserter( this->dim_lengths ));
	state = ALLOCATED;
}

//void CVariable_Container::alloc(int dims,  ...){
//	int i;
//	std::vector<int> lengths(dims);
//	std::va_list args;
//	va_start(args, dims);
//	for (i = 0; i < dims; i++) lengths[i] = va_arg(args, int);
//	va_end(args);
//	alloc(lengths);
//}

void CVariable_Container::add_existing_scalar(void* var){
	//adds a EXISTING scalar to the list of parameters.
	int size_elem = data_type.get_bytes();
	this->data = var;
	this->data_base = var;
	this->dims = 1;
	this->dim_lengths.clear();
	this->dim_lengths.push_back(1);
	this->existing = true; //this is the distinctive marker.
	//this->allocated = true;
	this->size_bytes = size_elem;
	this->size_elems = 1;
	state = ALLOCATED;
}

void CVariable_Container::allocate_space(const std::vector<int>& lengths){
	//validation
	std::string err_mess = name;
	switch (state) {
	case UNINITIALIZED:
		err_mess += " uninitialized";
		throw std::runtime_error(err_mess);
		break;
	case ALLOCATED:
		err_mess += " already allocated";
		throw std::runtime_error(err_mess);
		break;
	default:
		break;
	}
	if (lengths.size() != (unsigned int)this->dims){
		err_mess += ": dimensions do not match";
		throw std::runtime_error(err_mess);
	}
	//Everything clear. Just allocate.
	alloc(lengths);
}
void CVariable_Container::allocate_space(int _dims,...){
	std::vector<int> lengths(_dims);
	std::va_list args;
	va_start(args, _dims);
	for (int i = 0; i < _dims; i++) lengths[i] = va_arg(args, int);
	va_end(args);
	allocate_space(lengths);
}

void CVariable_Container::register_data(const std::vector<int>& dimensions, void* raw_data, bool cleanup){
	//validation
	std::string err_mess = name;
	switch (state) {
	case UNINITIALIZED:
		err_mess += ": uninitialized";
		throw std::runtime_error(err_mess);
		break;
	case ALLOCATED:
		err_mess += ": already allocated";
		throw std::runtime_error(err_mess);
		break;
	default:
		break;
	}
	if (dimensions.size() != (unsigned int)this->dims){
		err_mess += ": dimensions do not match";
		throw std::runtime_error(err_mess);
	}
	this->data = dan_flat2arrayND_cpp(raw_data, data_type.get_bytes(), dimensions);//if fails, throws bad_alloc
	this->data_base = raw_data;
	int size_elem = data_type.get_bytes();
	this->size_bytes = data_type.get_bytes(); 
	for (unsigned int i = 0; i < dimensions.size(); i++) size_bytes *= dimensions[i];
	this->existing = !cleanup; // <- controls if the object has to clean up after or not.
	this->size_elems = size_bytes / size_elem;
	std::copy(dimensions.begin(), dimensions.end(), std::back_inserter( this->dim_lengths ));
	state = ALLOCATED;
}

void CVariable_Container::swap_internal2external_scalar(void* existing_scalar){
	//Convert an internal scalar to an external storage
	std::string err_mess = this->name;
	if (state != ALLOCATED){
		err_mess += " not allocated";
		throw std::runtime_error(err_mess);
	} else if( this->dims != 1){
		err_mess += " not scalar";
		throw std::runtime_error(err_mess);
	} else if(this->existing){
		err_mess += " is externally allocated";
		throw std::runtime_error(err_mess);
	}
	//copy the current value of the internal storage to the external.
 	std::copy((char*)this->data_base, (char*)this->data_base + size_bytes, (char*)existing_scalar);
	// return the memory to the system.
	::operator delete (data_base); 
	add_existing_scalar(existing_scalar);
}

void CVariable_Container::copy_from_raw(void* src){
	abort_if_not_allocated();
	std::copy((char*)src, (char*)src + this->size_bytes, (char*)this->data_base);
}
void CVariable_Container::copy_raw_data(void* dst){
	abort_if_not_allocated();
	std::copy((char*)data_base, (char*)data_base + size_bytes, (char*)dst);
}

//Accessor that have to check the state of the object.
inline
void CVariable_Container::abort_if_not_allocated(){
	if (state != ALLOCATED) {
		std::string err = name + " not allocated";
		throw std::runtime_error(err);
	}
}

std::vector<int>& CVariable_Container::get_dim_lengths(){ 
	abort_if_not_allocated();
	return dim_lengths; 
}
void* CVariable_Container::get_data()	{ 
	abort_if_not_allocated();
	return data; 
}

void* CVariable_Container::get_data_base() { 
	abort_if_not_allocated();
	return data_base; 
}

int	CVariable_Container::get_size_bytes() { 
	abort_if_not_allocated();
	return size_bytes; 
}

int	CVariable_Container::get_size_elems() { 
	abort_if_not_allocated();
	return size_elems; 
}

bool CVariable_Container::get_existing()	{ 
	abort_if_not_allocated();
	return existing; 
}
