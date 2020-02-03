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

#ifndef _CPARAMS_GENERIC_H
#define _CPARAMS_GENERIC_H

//TO Do: Create functions for extracting and loading metadata for serialization.

#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <stdexcept>
#include "CVariable_Container.h"

typedef  CVariable_Container CPar_defs;

class CParams_generic {
protected:

public:
	typedef std::map<std::string, CPar_defs*> t_dict_type;
	typedef CParams_generic::t_dict_type::iterator iterator_dictionary_t;
	CParams_generic(){}
	virtual ~CParams_generic(){
		for (std::map<std::string, CPar_defs*>::iterator it = dict.begin(); it != dict.end(); it++)
			delete it->second;
	}
	CParams_generic& operator = (CParams_generic& orig){
		//caution!!! this doesn't check that the two objects are compatible. use with care.
		int blobsize = get_blob_size();
		unsigned char *buffer = new unsigned char[blobsize];
		orig.serialize_out(buffer);
		this->serialize_in(buffer);
		delete[] buffer;
		return(*this);
	}
	CPar_defs& operator [] (const std::string& key){
		return *dict[key];
	}
	CParams_generic(CParams_generic& orig){
		//This copy constructor doesn't take into consideration the "existing" flag in CPar_defs objects. It allocates
		// all variables.
		for (t_dict_type::iterator it = dict.begin(); it != dict.end(); it++){
			//Create each data object.
			CPar_defs& par = *(it->second);
			CPar_defs new_par(par);
			dict[par.get_name()] = &par;
		}
		//populate it.
		(*this) = orig; //use the asignment operator.
	}
	//functions for managing the collection
	void* add(const std::string &key, CPar_Data_Type::data_type_t type, const std::vector<int>& lengths);
	void* add(const std::string &key, CPar_Data_Type::data_type_t type, int dims, ...);
	void* add(CPar_defs& Existing_Variable);
	void* add_existing_scalar(const std::string &key, CPar_Data_Type::data_type_t type, void* var);
	void erase_variable(const std::string& key){
		if (!check_key(key)) return;
		CPar_defs* p = dict[key];
		delete p;
		dict.erase(key);
	}
	//functions for delayed allocation
	void add_no_alloc(const std::string& key, CPar_Data_Type::data_type_t type, int dims){
		CPar_defs* p = new CPar_defs(type, dims, key);
		dict[key] = p; 
	}
	void* alloc_delayed(const std::string& key, const std::vector<int>& lengths){
		CPar_defs* p = dict[key];
		p->allocate_space(lengths);
		return( p->get_data() );
	}
	//information about the object
	int get_nVars(){ return(dict.size()); }
	int get_blob_size();
	std::vector<std::string> get_registered_keys(){
		std::vector<std::string> v;
		for (t_dict_type::iterator it = dict.begin(); it != dict.end(); it++){
			v.push_back(std::string(it->first));
		}
		return(v);
	}
	//serialization functions
	void serialize_out(unsigned char *buffer);
	void serialize_in(unsigned char* buffer);
	bool check_key(const std::string &key){		return(dict.count(key) > 0);}
	CPar_defs& get_dataObject(const std::string &key){ return(*dict[key]);}
	std::map<std::string, CPar_defs*>& get_container(){	return(dict);	}	
private:
	t_dict_type dict;
};
#endif
