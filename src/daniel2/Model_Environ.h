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

#ifndef _MODEL_ENVIRON_H
#define _MODEL_ENVIRON_H

#include "CData.h"
#include "CParam.h"
#include "CChain.h"
#include "dan_sys.h"
#include <map>
#include <vector>
#include <string>
#include <algorithm> //For std::copy

class CMCMC_Trace {
	//terms:
	//	capacity		: number of samples that the tracer can store
	//	size			: number of samples currently stored in tracer
	//	capacity_elems	: capacity in number of elements (e.g. doubles) that the tracer can store for a key
	//	capacity_bytes	: same, but measured in bytes
	//	size_elems		: current number of elements (e.g. doubles) for a parameter currently stored in the tracer.
	//	size_bytes		: same in bytes.

public:
	CMCMC_Trace(CParam* par, const int& len) : par_container(par->get_container()), next_sample(0), len(len){
	}
	void set_trace(const std::string& key){
		if (!par_container.check_key(key)){
			//No such a parameter.
			std::string s = "Parameter " + key + " not found.";
			throw std::runtime_error(s);
		} 
		if (Trace_Container.check_key(key)){
			//Trace already created.
			std::string s = "Tracer for " + key + " already exists.";
			throw std::runtime_error(s);
		}
		std::vector<int> ln(par_container[key].get_dims() + 1);
		ln[0] = len;
		std::copy(par_container[key].get_dim_lengths().begin(), par_container[key].get_dim_lengths().end(), ln.begin() + 1);
		try{
			Trace_Container.add(key, par_container[key].get_data_type().get_data_type(), ln);
		} catch (const std::exception& e){
			std::string s = "Can't allocate trace for " + key + " : " + e.what();
			throw std::runtime_error(s);
		}
	}
	void set_trace(const std::string& key, void* src, int dims, int* lengths, CPar_Data_Type::data_type_t type){
		//This one is tricky. It sets a new tracer for a parameter that does not come from the par_container object.
		std::vector<int> ln(dims + 1);
		ln[0] = len;
		std::copy(lengths, lengths + dims, ln.begin() + 1);
		try{
			Trace_Container.add(key, type, ln);
		} catch (const std::exception& e){
			std::string s = "Can't allocate trace for " + key + " : " + e.what();
			throw std::runtime_error(s);
		}
		extra_sources[key] = src;
	}
	bool check_trace(const std::string& key){
		return(Trace_Container.check_key(key));
	}
	void _Update(){
		//This is a public method, but should be used with care. It's better to have it 
		//  only accessed in a controlled way.
		if (next_sample == len) return;
		std::map<std::string, CPar_defs*> &thismap = Trace_Container.get_container();
		for (std::map<std::string, CPar_defs*>::iterator it = thismap.begin(); it != thismap.end(); it++){
			//get &(dest[index][0][0]...[0]). this is quite tricky.
			int size = it->second->get_size_bytes() / len; //size of each sample in bytes
			const std::string &key = it->first;
			void *dst = (void *)((char*)it->second->get_data_base() + next_sample * size);
			void *src = NULL;
			if (par_container.check_key(key)){
				src = par_container[key].get_data_base();
			} else if(extra_sources.count(key) >0) {
				src = extra_sources[key];
			} else {
				//THIS WOULD BE AN INCONSISTENT STATE. MAYBE WANT TO THROW AN ERROR?
				continue;
			}

			std::copy((char*)src, (char*)src + size, (char*)dst);
		}
		++next_sample;
	}
	void Clear(){
		next_sample = 0;
	}
	void Copy_trace(const std::string& key, void *dst){
		//Careful: Can't use the container copy functions because
		// we have to only copy (size) elements.
		int _len = this->get_size_bytes(key);
		char* src = (char*)Trace_Container[key].get_data_base();
		std::copy(src, src + _len, (char*)dst);
	}
	int get_capacity() {
		return(len);
	}
	int get_size() {
		return( next_sample );
	}
	int get_size_bytes(const std::string& key){
		//container.get_get_size_bytes() / len = number of bytes of the data element.
		return Trace_Container[key].get_size_bytes() / len * next_sample;
	}
	int get_size_elems(const std::string& key){
		return Trace_Container[key].get_size_elems() / len * next_sample ;
	}
	int get_capacity_bytes(const std::string& key){
		return Trace_Container[key].get_size_bytes();
	}
	int get_capacity_elems(const std::string& key){
		return Trace_Container[key].get_size_elems();
	}
	void* get_raw(const std::string& key){
		return Trace_Container[key].get_data_base();
	}
	void* get_array(const std::string& key){
		return Trace_Container[key].get_data();
	}
	std::vector<int>& get_dims(const std::string& key){
		return Trace_Container[key].get_dim_lengths();
	}
	std::vector<std::string> get_trace_keys(){
		//return by value. A lot of overhead, but it's safer this way.
		const std::vector<std::string>& r = Trace_Container.get_registered_keys();
		return(r);
	}
	int get_elem_size (const std::string& key){
		return Trace_Container[key].data_type.get_bytes();
	}
	CPar_Data_Type& get_data_type(const std::string& key){
		return Trace_Container[key].get_data_type();
	}


protected:
	CParams_generic Trace_Container; //desctructors called automatically when the object is destroyed.
	CParams_generic& par_container;
	std::map<std::string, void*> extra_sources;
	int next_sample;
	const int len;
private:
	//disable complicated constructors.
	CMCMC_Trace(); //default constructor
  	CMCMC_Trace(const CMCMC_Trace&); //copy constructor
};

//Note that this doesn't depend on any R functionality.
// Keep it that way.
class CModel_Environ_Simple_base{
public:
	//WARNING: model_base HAS TO SET model_base after using register_model . The idea is to pass this
	// before initializing 
	typedef enum {OBJECT_CREATED, MODEL_REGISTERED, MODEL_INITIALIZED, UPDATING} object_state_t;
	//<<---- 2016-02-11: removed superfluous intialization of model_base. CHECK FOR SIDE EFECTS.
	CModel_Environ_Simple_base(CData* _data_base, CParam *_par_base, int len_buffer, int subsamp = 1, bool _del_objects = true)
			:	updating_output(true), state(OBJECT_CREATED), data_base(_data_base), par_base(_par_base),//, model_base(model_base), 
			tracing(false), thining(subsamp), del_objects(_del_objects), iteration(0), count_subsamp(0){
		try{
			trace = new CMCMC_Trace(par_base, len_buffer);
		} catch (const std::exception& e){
			std::string err = "Trace creation failed: ";
			err += e.what();
			throw std::runtime_error(err);
		}
	
	}
	void register_model(CChain* model_){ 
		model_base = model_;
		state = MODEL_REGISTERED;
	}
	//override these methods to add new parameters to trace.
	virtual ~CModel_Environ_Simple_base(){
		if (del_objects){
			delete data_base;
			delete par_base;
		}
		delete model_base;
		delete trace;
	}
	virtual void init_model(){	
		if (state != MODEL_REGISTERED) throw std::runtime_error("Cannot initialize now");
		model_base->Initializes();
		this->iteration = 0;
		this->count_subsamp = 0;
		this->trace->Clear();
		state = MODEL_INITIALIZED;
	}
	virtual std::vector<std::string> get_param_keys(){ 
		return par_base->get_var_list();
	}
	virtual void set_trace(const std::string& par_key){
		const std::vector<std::string>& lista = get_param_keys();
#ifdef __sun
		//Solaris has some quirks around the function std::count.
		if(std::find(lista.begin(), lista.end(), par_key) == lista.end() ) return;
#else
		if (std::count(lista.begin(), lista.end(), par_key) == 0) return;
#endif
		trace->set_trace(par_key);
	}

	///****                                              *****///
	///**** DO NOT OVERRIDE THE FOLLOWING PUBLIC METHODS *****///
	///****												 *****///
	
	//1) ------------------
	//	Actual enviornment level actions and stuff
	object_state_t get_model_state(){
		return(this->state);
	}
	void activate_updating_output(){
		updating_output = true;
	}
	void deactivate_updating_output(){
		updating_output = false;
	}
	bool get_updating_ouptut(){
		return(updating_output);
	}
	void activate_tracing(){
		tracing = true;
	}
	void deactivate_tracing(){
		tracing = false;
	}
	bool get_tracing(){
		return(tracing);
	}
	int get_thinning(){
		return(this->thining);
	}
	int get_iteration(){
		return(iteration);
	}
	void update_model(int iterations = 1){
		if (state != MODEL_INITIALIZED && state != UPDATING)
			throw std::runtime_error("Model not initialized");
		state = UPDATING;
		for (int i = 0; i < iterations; i++, iteration++){
			if (updating_output){
				DAN_PRINTF("\riteration = %d       ", this->iteration + 1);
			}
			model_base->Update();
			Other_Updating_Options();
			if (tracing && ++count_subsamp == thining){ 
				count_subsamp = 0;
				Other_Tracing_Options();
				trace->_Update();
			}
		}
	}
	void reset_subsamp(const int& new_ss){
		this->thining = new_ss;
		this->trace->Clear(); //reset_traces();
	}
	void recreate_trace(const int& len_buffer){
		CMCMC_Trace* new_trace;
		//this throws an exception when it fails. Catch at some place
		new_trace = new CMCMC_Trace(par_base, len_buffer);
		
		//copy trace parameters
		const std::vector<std::string>& keys = trace->get_trace_keys();//get_trace_keys();
		for (std::vector<std::string>::const_iterator it = keys.begin(); it != keys.end(); it++){
			new_trace->set_trace(*it);
		}
		delete trace;
		trace = new_trace;
	}

	//2) ------------------------------ 
	// Interface to parameters object
	CParams_generic& get_params(){
		return par_base->get_container();
	}
	CPar_defs& get_param_container(const std::string& par_key){
		return(par_base->get_container().get_dataObject(par_key));
	}

	bool check_param_key(const std::string& key){
		return (*par_base).check_key(key);
	}
	void* get_param(const std::string& key_data){ 
		return (*par_base)[key_data].get_data_base();
	}
	int get_blob_size(){
		return(par_base->get_blob_size());
	}
	void get_blob(void* buffer){
		par_base->get_blob((unsigned char*)buffer);
	}
	void set_blob(void* buffer){
		par_base->set_blob((unsigned char*)buffer);
	}
	void copy_param(const std::string& key, void* dst){
		(*par_base)[key].copy_raw_data(dst);
	}
	void write_param(const std::string& key, void* src){
		(*par_base)[key].copy_from_raw(src);
	}
	int get_param_size(const std::string& key){
		return (*par_base)[key].get_size_bytes();
	}
	//3)---------------------
	//	Interface to trace object
	CMCMC_Trace& Tracer (){ return *trace; }

	//4)---------------------
	// Interface to CChain object
	CChain& get_CChain(){return *model_base;}

protected:
	bool updating_output;
	object_state_t state;
	CData* data_base;
	CParam * par_base;
	CChain* model_base;
	CMCMC_Trace* trace;
	bool tracing;
	int thining;
	const bool del_objects;
	int iteration, count_subsamp;
	///----------------------------------------------------------------------------------///
	///*** USE THESE METHODS TO ADD FUNCTIONALITY TO THE CLASS THROUGH INHERITANCE  *****///
	///----------------------------------------------------------------------------------///
	virtual void Other_Updating_Options(){};
	virtual void Other_Tracing_Options(){};
};


//Again. This doesn't depend on any R or other functionality. Keep it that way.
template<typename T_Data, typename T_Model, typename T_Param>
class CModel_Environ_Simple : public CModel_Environ_Simple_base{
	//This template specializes Model_Environ_Simple_base to handle actual the actual classes
	// derived from CData, CParam and CChain. It links the child classes to the (parent) pointers in
	// CModel_Environ_Simple so that that class can use the parents' interface to do stuff.
public:
	CModel_Environ_Simple(T_Data* data_, T_Param* par_, int len_buffer, int subsamp = 1, 
		bool del_objects = true) 
			:	CModel_Environ_Simple_base(data_, par_, len_buffer, subsamp, del_objects), data(data_), par(par_)
	{
		try {
			model = new T_Model(data_, par_);
		} catch (const std::exception& e){
			std::string err = "CModel_Environ_Simple couldn't create Model Object: ";
			err += e.what();
			throw std::runtime_error(err);
		}
		CModel_Environ_Simple_base::register_model(model);
	}
	T_Data& get_data(){ return(*data); } 
	T_Param& get_params(){ return(*par); }
	T_Model& get_model(){ return(*model);}
protected:
	T_Data* data;
	T_Param* par;
	T_Model* model;
};

//--------------------------------------------------------------------------------------//
// The easiest way to add other functionals to a model environment is to derive a class
//  from CModel_Environ_Simple that extend the functionality by computing the required 
//  statistics and tracing them.
// - To trace those statistics one need to register them into this->trace using
//   trace->set_trace(<KEY>, <PTR_TO_DATA_SOURCE>, <n_DIMS>, <LENGTHS_OF_DIMS>, <TYPE>)
// - To compute the statistics, override the virtual functions Other_Tracing_Options() 
//   and Other_Update_Options(). These functions are called each tracing checkpoint 
//   and iteration, respectively.
//--------------------------------------------------------------------------------------//
//Example:
//template<class D, class M, class P>
//class miclass : public CModel_Environ_Simple<D, M, P>{
//public:
//	typedef CModel_Environ_Simple<D, M, P> mit_t;
//	typedef CDisclosure<P> midisc_t;
//	miclass(D* data, P* par, int popsize, int len_buffer, int BI=0, int subsamp=1, bool del_objects=true) 
//			:	mit_t(data, par, len_buffer, BI, subsamp, del_objects), 
//				disc(new midisc_t(par) ){
//		disc->Set_PopSize(popsize);
//		disc->load_data(data);
//		int ln[] = {1};
//		trace->set_trace("tau1", &tau1, 1, ln, CParams_generic::T_INT);
//		trace->set_trace("tau2", &tau2, 1, ln, CParams_generic::T_DOUBLE);
//	}
//	~miclass(){
//		delete disc;
//	}
//protected:
//	int tau1;
//	double tau2;
//	void Other_Tracing_Options(){
//		disc->compute_prob_uniques();
//		disc->sample_synth_pop(model->r);
//		tau1 = disc->get_tau_sim();
//		tau2 = disc->get_tau2_sim();
//	}
//	void Other_Updating_Options(){
//		cout << "Iter = " << this->iteration <<  "tau1 = " << tau1 << endl;
//	}
//	midisc_t *disc;
//};

#endif
