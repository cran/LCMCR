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

#include <cstdlib>
#include <cstdarg> 
#include "dan_array_utils.h"
#include "CParams_generic.h"
#include <iostream>
#include <iterator>
#include <string>
#include <stdexcept>




///////////////////////////////////////////////////
//Implementation of CParams_generic
///////////////////////////////////////////////////

void* CParams_generic::add(CPar_defs& Existing_Variable){
	if (check_key(Existing_Variable.get_name()) ){
		throw std::runtime_error("Variable already exists");
	}
	dict[ Existing_Variable.get_name() ] = &Existing_Variable;
	return Existing_Variable.get_data();
}


void* CParams_generic::add(const std::string &key, CPar_Data_Type::data_type_t type, const std::vector<int>& lengths){
	//allocates space and adds new key.
	if (lengths.size()< 1) return NULL;
	//size_elem = type_size(type);
	CPar_defs* descr;
	try{
		descr = new CPar_defs(type, lengths, key);
	} catch (...){
		std::string m = "Can't create parameter object \"" + key + "\"";
		throw std::runtime_error(m);
	}
	#pragma omp critical
	dict[key] = descr;
	return descr->get_data();
}

void* CParams_generic::add(const std::string &key, CPar_Data_Type::data_type_t type, int dims,  ...){
	int i;
	std::vector<int> lengths(dims);
	std::va_list args;
	va_start(args, dims);
	for (i = 0; i < dims; i++){
		lengths[i] = va_arg(args, int);
	}
	va_end(args);
	return add(key, type, lengths);
}

void* CParams_generic::add_existing_scalar( const std::string &key, CPar_Data_Type::data_type_t type, void* var){
	//adds a EXISTING scalar to the list of parameters.
	//int size_elem = type_size(type);//(type==T_INT) ? sizeof(int) : sizeof(double);
	CPar_defs* descr = new CPar_defs(type, var, key);
	dict[key] = descr;
	return var;
}

//void* CParams_generic::get(const std::string &key){
//	return dict[key]->get_data();
//}

//void* CParams_generic::get_base_address(const std::string &key){
//	return dict[key]->get_data_base();
//}

//int CParams_generic::get_size_bytes(const std::string &key){
//	return dict[key]->get_size_bytes();
//}

int CParams_generic::get_blob_size(){
	int size = 0;
	std::map<std::string, CPar_defs*>::iterator it;
	for (it = dict.begin(); it != dict.end(); it++){
		size += it->second->get_size_bytes();
	}
	return(size);
}

void CParams_generic::serialize_out(unsigned char *buffer){
	typedef unsigned char* pt_byte;		
	pt_byte ptr_buff = buffer;
	std::map<std::string, CPar_defs*>::iterator it;
	for (it = dict.begin(); it != dict.end(); it++){
		std::copy(pt_byte(it->second->get_data_base()), pt_byte(it->second->get_data_base()) + it->second->get_size_bytes(), ptr_buff);
		ptr_buff += it->second->get_size_bytes();
	}
}

void CParams_generic::serialize_in(unsigned char *buffer){
	typedef unsigned char* pt_byte;		
	pt_byte ptr_buff = buffer;
	std::map<std::string, CPar_defs*>::iterator it;
	for (it = dict.begin(); it != dict.end(); it++){
		std::copy(ptr_buff, ptr_buff + it->second->get_size_bytes(), pt_byte(it->second->get_data_base()));
		ptr_buff += it->second->get_size_bytes();
	}
}


//////////
