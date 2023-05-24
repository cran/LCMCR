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

#ifndef _CLcm_Bay_BASIC_FREQ_H
#define _CLcm_Bay_BASIC_FREQ_H

#include <stdio.h>
#include <stdlib.h>
#include "daniel2/CChain.h"
#include "daniel2/CParam.h"
#include "CData_DM.h"


class CParams_NPLCM_CR_Basic_Freq : public CParam {
public:
	//parameters
	int J, K, n, M; //C=levels of covariate (number of strata)
	//int *zI; //(z_i = k) i=1.. n
	int **count_zIK;
	//double **lambdaJK; //(psi_{jkl}) l -> # of levels
	double ***log_lambdaJK2; //(psi_{jkl}) l -> # of levels
	double *nuK, *log_nuK;
	int *countK;
	int *count0K;  //for unobserved individuals
	double alpha;
	//augmented
	int n0; 
	//int *z0I;
	double prob_zero;
	//auxiliary
	int*** aux_JK2; //to store #{i: cov_i=c, z_i=k, x_{ij}=x}
	int k_star;
	//priors
	double a_alpha, b_alpha; // alpha ~ Gamma[a_alpha, b_alpha]

	//Constructors.
	CParams_NPLCM_CR_Basic_Freq(int _J, int _K, int _n, int _M, double _a_alpha, double _b_alpha)
		:	J(_J), K(_K), n(_n), M(_M), a_alpha(_a_alpha), b_alpha(_b_alpha) {
		class_construct();
	}
	CParams_NPLCM_CR_Basic_Freq(CData_DM* _dat, int _K, double _a_alpha, double _b_alpha)
		:	J(_dat->J), K(_K), n(_dat->n), M(_dat->ncells), a_alpha(_a_alpha), b_alpha(_b_alpha) {
		class_construct();		
	}
	CParams_NPLCM_CR_Basic_Freq(CParams_NPLCM_CR_Basic_Freq& orig)
		:	J(orig.J), K(orig.K), n(orig.n), M(orig.M), a_alpha(orig.a_alpha), b_alpha(orig.b_alpha){
		class_construct();
		//*this = orig; //copy the container. 
		this->storage = orig.storage; 
	}
	CParams_NPLCM_CR_Basic_Freq& operator=(CParams_NPLCM_CR_Basic_Freq& orig){
		CParam::operator=(orig);
		return(*this);
	}

	void initizalize(gsl_rng *r); //Initialize the parameters
private:
	void class_construct();
};
/////////////////////////////////
/////////////////////////////////

class CNPLCM_CR_Basic_Freq : public CChain{
public:
	//CNPLCM_CR(){}; //default constructor. DO NOT USE!!!!
	CNPLCM_CR_Basic_Freq(CData_DM *_data, CParams_NPLCM_CR_Basic_Freq*_par) : CChain(_par), par(_par) , data(_data){
		class_construct(data, par);
	}
	CNPLCM_CR_Basic_Freq(CData_DM* _data, int _K, double _a_alpha, double _b_alpha) 
			:	CChain(), data(_data){
		
				register_param(par = new CParams_NPLCM_CR_Basic_Freq(data->J, _K, data->n, data->ncells, _a_alpha, _b_alpha));
		class_construct(data, par);
	}
	void Update();
	CParams_NPLCM_CR_Basic_Freq *par; 
	void WriteInfoSheet(char *filename);
	void Initializes();
	void Get_Probs_Observed_Data(double *probs_n);
private:
	CData_DM* data;
	void class_construct(CData_DM* d, CParams_NPLCM_CR_Basic_Freq* p){
		this->K = p->K; J = d->J; n = d->n; 
		setModelSignature("CNPLCM_CR_Basic_Freq");
	}
protected:
	int J, n; //to be set with values from *data.
	int K; //to be set with values from *par or in constructor.
	//Full conditional samplers
	//  all these functions update the values on the *par object.
	void sam_nu();
	void sam_lambda();
	void sam_countzIK();
	void sam_alpha();
	
	//augmented sample:
	void sam_n0();
	void sam_z0x0();

	//auxiliary functions
	void permute_latent_classes_by_nu();
	void compute_probs_miss();
	void tabulate_elements();
	void CountKs();
};

#endif  //_CLcm_Bay_CR_H
