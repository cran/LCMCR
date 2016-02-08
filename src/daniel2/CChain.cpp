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


#include "CChain.h"
#include "dan_math_gsl.h"
#include <string.h>
#include <time.h>

//Constructor
void CChain::class_construct(){
	this->verbose = true;
	this->r = gsl_rng_alloc(gsl_rng_taus2);
	int off = sizeof(time_t) / 2;
	int t = ( int(time(0)) << off) >> off;
	gsl_rng_set(r,t);
	current_iteration = 0;
	this->start_time = clock(); //just to set to something. Modify it with "reseed_rng"
	setModelSignature("MODEL_SIGNATURE_NOT_SET"); //HAS TO BE SET IN CHILD CLASS!!!!!!!
}

void CChain::Initializes(){
	param_parent->initizalize(r);
}

void CChain::Update(){
	this->current_iteration++;
	lap_clock(); //Reset the timer
	//Add functionality in the child class. Don't forget to run this part of the code.
}

void CChain::lap_clock(){
	long stop_time;
	stop_time = clock();
	this->lap_time = double(stop_time-start_time)/double(CLOCKS_PER_SEC);
	start_time = stop_time;
}

void CChain::setModelSignature(const std::string &signature){
	model_signature.assign(signature);
}

const std::string&  CChain::getModelSignature(){
	return(this->model_signature);
}
