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

/*
	CChain class definition
	Generic interface class for a MCMC object. Methods have to be implemented in 
	 classes inheriting from this one.

*/

#if !defined(_CCHAIN_H)
#define _CCHAIN_H

#include <gsl/gsl_rng.h>
#include "CParams_generic.h"
#include "CParam.h"
#include <string.h>
#include <string>
#include <time.h>
#include <string>

/*
	Child classes have to initialize param_parent to an instance of a class inherited from CParams.
*/

class CChain {
	//1) Child class must create and register an object derived from CParam: 
	//		either with register_param(par = new CPara...(a,b,c...)); or directly using the second constructor.
	//2) Child class must set model signature using setModelSignature(const std::string &);
	//3) Do not allocate any structure that defines the state of the model here. Use the CParam object.
	//4) Do not use the members commented as "deprecated" for future projects.
public:
	//Constructors. 
	CChain() { class_construct();}
	CChain(CParam *par) : param_parent(par) {class_construct();} 
	virtual ~CChain(){};
	void setModelSignature(const std::string &signature);
	const std::string& getModelSignature();
	virtual void Initializes() = 0;
	virtual void Update()=0;
	double get_lap_time(){ return(lap_time);}
	virtual void reseed_rng(const unsigned int new_seed){ gsl_rng_set(r, new_seed); }
	void messages_off(){ verbose = false;}
	void messages_on(){ verbose = true;}
	gsl_rng *r;

protected:
	CParam* param_parent; //may need to be public
	int unsigned long ran_seed;
	double lap_time;// time between two consecutive calls to lap_clock();
	void register_param(CParam* p){param_parent = p;}
	void lap_clock();
	int current_iteration;
	bool verbose; //Write messages to the console?
	CParams_generic local_vars;//Local allocation ***THE STATE OF THE CHAIN MUS NOT DEPEND OF THESE VALUES***
	void* add_static(const std::string& key,
						CPar_Data_Type::data_type_t type,
						const std::vector<int>& lengths) {
		return(local_vars.add_static(key, type, lengths));
	}
	void* add_static(const std::string& key,
						CPar_Data_Type::data_type_t type, int dims, ...) {
							//If key exists, just return a pointer to that data
		std::vector<int> lengths(dims);
		std::va_list args;
		va_start(args, dims);
		for (int i = 0; i < dims; i++) {
		lengths[i] = va_arg(args, int);
		}
		va_end(args);
		return local_vars.add_static(key, type, lengths);
	}
private:
	std::string model_signature; //MUST BE INTIALIZED IN CHILD CLASS.
	long start_time;
	void class_construct();
};

#endif  //_CCHAIN_H
