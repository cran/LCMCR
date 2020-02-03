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

//--------------------------------------------------------------------------------
// LCMCR with incomplete stratification
//--------------------------------------------------------------------------------
#ifndef _CLcm_Bay_STRAT_FREQ_H
#define _CLcm_Bay_STRAT_FREQ_H

#include <stdio.h>
#include <stdlib.h>
#include "daniel2/CChain.h"
#include "daniel2/CParam.h"
#include "CData_DM_Strat.h"


class CParams_NPLCM_CR_Strat_Freq : public CParam {
public:
	//Hyperpriors
	int K;
	double a_alpha, b_alpha; // alpha ~ Gamma[a_alpha, b_alpha]
	double a_lambda, b_lambda; // lambda_cjk ~ Beta(a_lambda, b_lambda)
	double N_Max_factor; // upper truncation multiplier (upper limit as a factor of observed).

	//For generation of labels
	int ***count_zmisCMK;
	int ***count_zstratCMK;
	int **count_z0CK;  //for unobserved individuals

	//For tabulation
	int **count_CK;
	int**** count_CJK2; //to store #{i: cov_i=c, z_i=k, x_{ij}=x}

	//For imputation of labels
	int **n_misMC; //imputed individuals with pattern m in stratum c.
	int *n_misC; //imputed individuals in stratum c.

	//For generation of counts
	double *prob_zeroC;
	int *n0C; 

	//Parameters
	//double ***lambdaCJK; 
	double ****log_lambdaCJK2; 
	double **nuCK, **log_nuCK;
	double *alphaC;
	double *rhoC;

	//Active lists:
	int min_n;
	int min_lists;
	int **activeCJ;
	int *activeC;
	//Constructors.

	CParams_NPLCM_CR_Strat_Freq(
		CData_DM_Strat *_dat,
		int _K, double _a_alpha, double _b_alpha, double _a_lambda, double _b_lambda,
		int _min_n, int _min_lists, double _N_Max_factor)
		: K(_K),
		  a_alpha(_a_alpha), b_alpha(_b_alpha),
		  a_lambda(_a_lambda), b_lambda(_b_lambda),
		  N_Max_factor(_N_Max_factor),
		  min_n(_min_n), min_lists(_min_lists)
	{
		class_construct(_dat);		
	}
	//CParams_NPLCM_CR_Strat_Freq(CParams_NPLCM_CR_Strat_Freq& orig)
	//	:	 K(orig.K), a_alpha(orig.a_alpha), b_alpha(orig.b_alpha){
	//	class_construct();
	//	this->storage = orig.storage; 
	//}
	CParams_NPLCM_CR_Strat_Freq& operator=(CParams_NPLCM_CR_Strat_Freq& orig){
		CParam::operator=(orig);
		return(*this);
	}

	void initizalize(gsl_rng *r); //Initialize the parameters
protected:
	void determine_active_lists(CData_DM_Strat* _dat, const int& min_n, const int& min_lists);
private:
	CParams_NPLCM_CR_Strat_Freq(CParams_NPLCM_CR_Strat_Freq& orig){}; //disable copy constructor.
	CParams_NPLCM_CR_Strat_Freq(){}; //disable default constructor
	void class_construct(CData_DM_Strat* dat);


};
/////////////////////////////////
/////////////////////////////////

class CNPLCM_CR_Strat_Freq : public CChain {
public:
	CParams_NPLCM_CR_Strat_Freq *par; 
	CNPLCM_CR_Strat_Freq(CData_DM_Strat *_data, CParams_NPLCM_CR_Strat_Freq*_par) : CChain(_par), par(_par) , data(_data){
		class_construct(data, par);
	}
	CNPLCM_CR_Strat_Freq(
		CData_DM_Strat *_data, int _K,
		double _a_alpha, double _b_alpha, double _a_lambda, double _b_lambda,
		int _min_n, int _min_lists, double _N_Max_factor)
		: CChain(), data(_data)
	{

		register_param(
			par = new CParams_NPLCM_CR_Strat_Freq(
				_data, _K, _a_alpha, _b_alpha, _a_lambda, _b_lambda,
				_min_n, _min_lists, _N_Max_factor)
		);
		class_construct(data, par);
	}
	void Update();
	void WriteInfoSheet(char *filename);
	void Initializes();
	//void Get_Probs_Observed_Data(double *probs_n);
private:
	CData_DM_Strat* data;
	CNPLCM_CR_Strat_Freq(){}; //Disable default constructor.
	void class_construct(CData_DM_Strat* d, CParams_NPLCM_CR_Strat_Freq* p){
		setModelSignature("CNPLCM_CR_Strat_Freq");
	}
	bool _is_zeros_strat(int* X, int c){
		bool ret = true;
		for (int j = 0; j < data->J - 1; j++){
			if (!par->activeCJ[c][j]) continue;
			if (X[j] == 1){
				ret = false;
				break;
			}
		}
		return(ret);
	}
protected:
	//Full conditional samplers
	//  all these functions update the values on the *par object.
	void sam_lambdaCJK();
	void sam_nuCK();
	void sam_alphaC();
	void sam_rhoC();

	//imputation of latent classes
	void sam_zmisCMK();
	void sam_zstratCMK();
	void sam_z0CK();

	//imputation of strata labels
	void sam_n_misMC();
	
	//augmented sample:
	void compute_probs_miss();
	void sam_n0C();

	//Aux function
	void reset_counters();
};

#endif  //Inclusion guard
