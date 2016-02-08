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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h> //special functions (for gamma)
#include <string.h>
#include <math.h>
#include "daniel2/dan_math_gsl.h"
#include "definitions.h"
#include "CData_DM.h"
#include "NPLCM_CR_Basic_Freq.h"
#include <vector>
#include <algorithm>
/*------------------------------------------------------------
Implementation of class CParams_NPLCM_CR 
------------------------------------------------------------*/
//Generic Constructor
void CParams_NPLCM_CR_Basic_Freq::class_construct(){
	//allocate
	_add_parameter("lambdaJK", CPar_Data_Type::T_DOUBLE, &lambdaJK, 2, J, K);
	//_add_parameter("zI",  CPar_Data_Type::T_INT, &zI, 1, n);
	_add_parameter("count_zIK",  CPar_Data_Type::T_INT, &count_zIK, 2, M, K);
	_add_parameter("nuK", CPar_Data_Type::T_DOUBLE, &nuK, 1, K);
	_add_parameter("log_nuK", CPar_Data_Type::T_DOUBLE, &log_nuK, 1, K);
	_add_parameter("countK", CPar_Data_Type::T_INT, &countK, 1, K);
	_add_parameter("count0K", CPar_Data_Type::T_INT, &count0K, 1, K);
	_add_parameter("aux_JK2", CPar_Data_Type::T_INT, &aux_JK2, 3, J, K, 2); 
	// register existing
	_add_existing_scalar("K", CPar_Data_Type::T_INT, &K);
	_add_existing_scalar("k_star", CPar_Data_Type::T_INT, &k_star);
	_add_existing_scalar("a_alpha", CPar_Data_Type::T_DOUBLE, &a_alpha);
	_add_existing_scalar("b_alpha", CPar_Data_Type::T_DOUBLE, &b_alpha);
	_add_existing_scalar("alpha", CPar_Data_Type::T_DOUBLE, &alpha);
	_add_existing_scalar("n0", CPar_Data_Type::T_INT, &n0);
//_add_existing_scalar("n0_max", CPar_Data_Type::T_INT, &n0_max);
	_add_existing_scalar("prob_zero", CPar_Data_Type::T_DOUBLE, &prob_zero);
	_add_existing_scalar("M", CPar_Data_Type::T_INT, &M);
}

void CParams_NPLCM_CR_Basic_Freq::initizalize(gsl_rng *r){
	//intialize parameters at reasonable values.
	alpha = 1.0; //initialize at full spread: uniform prior.
	prob_zero = 0.0;
	n0 = 0;
	k_star = 0;
	(*this)["nuK"].fill(1.0/double(K));
	(*this)["log_nuK"].fill(-log(double(K)));
	(*this)["lambdaJK"].fill(0.5);
	(*this)["countK"].fill(0);
	(*this)["count0K"].fill(0);
	(*this)["count_zIK"].fill(0); 
	//(*this)["zI"].fill(0);
	(*this)["aux_JK2"].fill(0);
}


//--------------------------------------------------------------------------------
// Implementation of class CNPLCM_CR
//--------------------------------------------------------------------------------
void CNPLCM_CR_Basic_Freq::Update() {
	//sam_z();
	sam_countzIK();
	sam_lambda(); 
	sam_nu();
	//imputation
	compute_probs_miss();
	sam_n0();
	sam_z0x0();
	//regularization
	sam_alpha();
}

void CNPLCM_CR_Basic_Freq::Get_Probs_Observed_Data(double *probs_n){
	//compute the probability of each observed data point.
	// user must provide properly allocated space.
	for (int i = 0; i < data->n; i++){
		double sum = 0.0;
		for (int k = 0; k < par->K; k++){
			double prod = par->nuK[k];
			for (int j = 0; j < par->J; j++){
				prod *=  par->lambdaJK[j][k];
			}
			sum += prod;
		}
		probs_n[i] = sum;
	}
}
 

inline
void CNPLCM_CR_Basic_Freq::sam_lambda(){
	int tmp1 = 1, tmp0 = 0;
	CVariable_Container& lJK = (*par)["lambdaJK"];
	CVariable_Container& raux_JK2 = (*par)["aux_JK2"];
	this->tabulate_elements(); //accumulate the counts of samples.
	//add the prior (+1) and sample
	for(int j = 0; j < data->J; j++){
		for(int k = 0; k < par->K; k++){
			//add prior and cast and sample!
			//par->lambdaJK[j][k] = gsl_ran_beta(r, double(1 + par->aux_JK2[j][k][1]), double(1 + par->aux_JK2[j][k][0]));
			lJK.arr_elem<double>(j,k) = gsl_ran_beta(r, double(1 + raux_JK2.arr_elem<int>(j, k, tmp1)), 
				double(1 + raux_JK2.arr_elem<int>(j, k, tmp0)));
		}
	}
}

//inline
//void CNPLCM_CR_Basic_Freq::sam_z(){
//	int i, j, k;
//	double prod;
//	static double probs[1000];
//	for (i = 0; i < data->n; i++){
//		for (k = 0; k < par->K; k++){
//			prod= par->nuK[k];
//			for (j = 0; j < par->J; j++){
//				prod *= data->x[i][j] == 1 ? par->lambdaJK[j][k] : 	1.0 - par->lambdaJK[j][k];
//			}
//			probs[k] = prod;
//		}
//		par->zI[i] = dan_multinomial_1(r, par->K, probs, false);
//	}
//}

inline
void CNPLCM_CR_Basic_Freq::sam_countzIK(){
	double probs[1000];//should have K elements
	double prod = 1.0;
	for (int m = 0; m < data->ncells; m++){
		for (int k = 0; k < par->K; k++){
			prod= par->nuK[k];
			for (int j = 0; j < par->J; j++){
				prod *= data->cells[m][j] == 1 ? par->lambdaJK[j][k] : 	1.0 - par->lambdaJK[j][k];
			}
			probs[k] = prod;
		}
		//sample!
		if (data->freqM[m] > 1){
			gsl_ran_multinomial(r, par->K, data->freqM[m], probs, (unsigned int*)(par->count_zIK[m]));
		} else {
			int z = dan_multinomial_1(r, par->K, probs, false);
			std::fill(par->count_zIK[m], par->count_zIK[m] + par->K, 0);
			par->count_zIK[m][z] = 1;
		}
	}
}

inline
void CNPLCM_CR_Basic_Freq::CountKs(){
	std::copy(par->count0K, par->count0K + par->K, par->countK);
	for (int m = 0; m < par->M; m++){
		for (int k = 0; k < par->K; k++){
			//par->countK[ par->zI[i] ]+ = ;
			par->countK[k] += par->count_zIK[m][k];
		}
	}
	par->k_star = par->K - std::count(par->countK, par->countK + par->K, 0);
}



inline
void CNPLCM_CR_Basic_Freq::sam_nu(){
	int k;
	double b = 0.0, a = 0.0, lgamma1, lgamma2, lsumgamma;
	int n_acc = 0;
	double l_acc_prod = 0.0;
	CountKs();
	for (k = 0; k < par->K - 1; ++k){
		n_acc += par->countK[k];
		a = double(1 + par->countK[k]);
		b = par->alpha + double(par->n0+ data->n - n_acc);
		lgamma1 = dan_log_gamma_1(r, a);
		lgamma2 = dan_log_gamma_1(r, b);
		lsumgamma = dan_log_sum(lgamma1, lgamma2);
		par->log_nuK[k] = (lgamma1 - lsumgamma + l_acc_prod);
		l_acc_prod += lgamma2 - lsumgamma; //log(1-V)
		par->nuK[k] = exp(par->log_nuK[k]);
	}
	par->log_nuK[par->K -1] = l_acc_prod;
	par->nuK[par->K - 1] = exp(par->log_nuK[par->K-1]);
}


inline
void CNPLCM_CR_Basic_Freq::sam_alpha(){
	par->alpha = gsl_ran_gamma(r, par->a_alpha + double(par->K) - 1.0, 1.0 / (par->b_alpha - par->log_nuK[par->K-1]));
}

inline
void CNPLCM_CR_Basic_Freq::compute_probs_miss(){
	int k, j;
	par->prob_zero = 0.0;
	for (k = 0; k < par->K; ++k){
		double prod = par->nuK[k];
		for (j = 0; j < data->J; ++j){
			prod *= 1.0 - par->lambdaJK[j][k];
		}
		par->prob_zero += prod;
	}
}


inline
void CNPLCM_CR_Basic_Freq::sam_n0(){
	//GSL: "This function returns a random integer from the negative binomial distribution, 
	// the number of failures occurring before n successes in independent trials with probability 
	// p of success. The probability distribution for negative binomial variates is,
	// p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k NOTE THE +1 IN GAMMA(K+1)!!!!
	par->n0 = gsl_ran_negative_binomial(r, 1.0 - par->prob_zero, data->n); //<-=-this one is the correct.
	//This corresponds to an (improper) prior P(Nmis|n) \propto n/(Nmis+n). This ensures that after integrating Nmis
	// the joint probability correspond to the correct truncated distribution.
}

inline
void CNPLCM_CR_Basic_Freq::sam_z0x0(){
	std::vector<double> table_z(par->K);
	std::vector<unsigned int> counts(par->K);
	for (int k = 0; k < par->K; k++){
		double p = par->nuK[k] ;
		for(int j = 0; j < par->J; j++){
			p *= 1.0 - par->lambdaJK[j][k];
		}
		table_z[k] = p;
	}
	gsl_ran_multinomial(r, par->K, par->n0, &(table_z[0]), (unsigned int*)(par->count0K));
}

inline
void CNPLCM_CR_Basic_Freq::tabulate_elements(){
	//this can be optimized by tabulating at the time of sampling z, z0, cov, and cov0.
	//CVariable_Container& aux_ref = (*par)["aux_JK2"];
	CVariable_Container& aux_ref = (*par)["aux_JK2"];
	aux_ref.fill(0);
	//std::fill(par->aux_JK2[0][0], par->aux_JK2[0][0] + par->J *par->K * 2, 0);
	//for (int i = 0; i < par->n; i++){
	//	int z = par->zI[i];
	//	int *xiJ = data->x[i];
	//	for (int j = 0; j <  par->J; j++){
	//		//par->aux_JK2[j][z][ xiJ[j] ]++;
	//		aux_ref.arr_elem<int>(j, z, xiJ[j])++;
	//	}
	//}
	for (int m = 0; m < par->M; m++){
		int *cell = data->cells[m];	
		for (int k = 0; k < par->K; k++){
			for (int j = 0; j <  par->J; j++){
				//par->aux_JK2[j][k][ cell[j] ] += data->freqM[m];
				aux_ref.arr_elem<int>(j, k, cell[j]) += par->count_zIK[m][k]; //data->freqM[m];
			}
		}
	}
	////augmented section
	for (int j = 0; j < par->J; j++){
		for (int k = 0; k < par->K; k++){
			//par->aux_JK2[j][k][0] += par->count0K[k];
			int tmp = 0;
			aux_ref.arr_elem<int>(j, k, tmp) += par->count0K[k];
		}
	}
}



void CNPLCM_CR_Basic_Freq::Initializes(){
	double g1, g2, gs;
	CChain::Initializes();

	//init nu and alpha
	par->alpha = 1;
	double lacc = 0.0;
	for (int k = 0; k < par->K - 1; k++){
		g1 = dan_log_gamma_1(r, 1);
		g2 = dan_log_gamma_1(r, 1);
		gs = dan_log_sum(g1, g2);
		par->log_nuK[k] = g1 - gs + lacc;
		lacc += g2 - gs;
	}
	par->log_nuK[par->K -1] = lacc;
	par->nuK[par->K-1] = exp(lacc);
	//init lambda with marginal relative frequencies (by stratum).
	std::vector<int> sum_J(par->J);
	std::fill(sum_J.begin(), sum_J.end(), 0);
	for (int j = 0; j < par->J; j++){
		//for (int i = 0; i < par->n; i++){
		//	sum_J[j] += data->x[i][j];
		//}
		for (int m = 0; m < par->M; m++){
			sum_J[j] += data->cells[m][j]  == 1 ? data->freqM[m] : 0;
		}		
		double p = double(sum_J[j]) / double(par->n);
		for (int k = 0; k < par->K; k++){
			par->lambdaJK[j][k] = p;
		}
	}
	par->alpha = 0.3 / double(K);
	for(int k =0; k < 4; k++){
		par->nuK[k] = 0.25 * 0.9;
	}
	for(int k = 4; k < par->K; k++){
		par->nuK[k] = 0.1 / double (par->K -4);
	}
	//warm-up run. In order to fix the 
	DAN_PRINTF("WARMING UP...\n");
	for (int i = 0; i < 500; ++i){
		//sam_z();
		sam_countzIK();
		sam_lambda();
		sam_nu();
		//truncation!
		compute_probs_miss();
		sam_n0();
		sam_z0x0();
		//and regularization!
		//sam_alpha();
		//DAN_PRINTF("\riter = %d kstar = %d \n", i, par->k_star);
	}
	//reorder latent classes
	permute_latent_classes_by_nu();
}

typedef std::pair<double, int> _CNPLCM_CR_Basic_Freq_mypair;
bool _CNPLCM_CR_Basic_Freq_comparator (const _CNPLCM_CR_Basic_Freq_mypair& l, const _CNPLCM_CR_Basic_Freq_mypair& r){
	return (l.first < r.first);
}
void CNPLCM_CR_Basic_Freq::permute_latent_classes_by_nu(){
	//reorder latent classes increasingly according to nuK. 
	std::vector<int> new2old(par->K), old2new(par->K);
	std::vector<_CNPLCM_CR_Basic_Freq_mypair> s_index(par->K);
	for (int k = 0; k < par->K; k++){
		s_index[k].first = 1.0 - par->nuK[k];
		s_index[k].second = k;
	}
	std::sort(s_index.begin(), s_index.end(), _CNPLCM_CR_Basic_Freq_comparator);
	std::vector<int> perm(K);
	for (int k = 0; k < par->K; k++){
		new2old[k] = s_index[k].second;
	} // new2old[a] == b <=> former b-th element is now a-th element.
	for (int k = 0; k < par->K; k++){
		old2new[ new2old[k] ] = k;
	}
	CParams_NPLCM_CR_Basic_Freq *tmp = new CParams_NPLCM_CR_Basic_Freq(*par);
	//for (int i = 0; i < par->n; i++){
	//	tmp->zI[i] = old2new[ par->zI[i] ];
	//}
	for (int m = 0; m < par->M; m++){
		for (int k = 0; k < par->K; k++){
			tmp->count_zIK[m][k] = par->count_zIK[m][ new2old[k] ];
		}
	}	
	for (int j = 0; j < par->J; j++){
		for (int k = 0; k < par->K; k++){
			tmp->lambdaJK[j][k] = par->lambdaJK[j][ new2old[k] ];
			tmp->aux_JK2[j][k][0] = par->aux_JK2[j][ new2old[k] ][0];
			tmp->aux_JK2[j][k][1] = par->aux_JK2[j][ new2old[k] ][1];
		}
	}
	for (int k = 0; k < par->K; k++){
		tmp->nuK[k] = par->nuK[new2old[k]];
		tmp->log_nuK[k] = par->log_nuK[new2old[k]];
		tmp->countK[k] = par->countK[new2old[k]];
		tmp->count0K[k] = par->count0K[new2old[k]];
	}
	//std::copy(tmp->zI, tmp->zI + par->n, par->zI);
	std::copy(tmp->count_zIK[0], tmp->count_zIK[0] + par->M * par->J, par->count_zIK[0]);
	std::copy(tmp->lambdaJK[0], tmp->lambdaJK[0] + par->J * par->K, par->lambdaJK[0]);
	std::copy(tmp->aux_JK2[0][0], tmp->aux_JK2[0][0] + par->J * par->K * 2, par->aux_JK2[0][0]);
	std::copy(tmp->nuK, tmp->nuK + par->K, par->nuK);
	std::copy(tmp->log_nuK, tmp->log_nuK + par->K, par->log_nuK);
	std::copy(tmp->countK, tmp->countK + par->K, par->countK);
	std::copy(tmp->count0K, tmp->count0K + par->K, par->count0K);
	delete tmp;
}

