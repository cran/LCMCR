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
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <string>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h> 
#include "daniel2/dan_math_gsl.h"
#include "CData_DM.h"
#include "NPLCM_CR_Basic_Freq.h"
#include "CData_DM_Strat.h"
#include "NPLCM_CR_Strat_Freq.h"

/*------------------------------------------------------------
Implementation of class CParams_NPLCM_CR 
------------------------------------------------------------*/
//Generic Constructor
void CParams_NPLCM_CR_Strat_Freq::class_construct(CData_DM_Strat* _dat){
	
	_add_parameter("activeCJ", CPar_Data_Type::T_INT, &activeCJ, 2, _dat->C, _dat->J - 1); 
	_add_parameter("activeC", CPar_Data_Type::T_INT, &activeC, 1, _dat->C); 

	//allocate the rest
	_add_parameter("n_misMC", CPar_Data_Type::T_INT, &n_misMC, 2, _dat->ncells_mis, _dat->C);
	_add_parameter("n_misC", CPar_Data_Type::T_INT, &n_misC, 1, _dat->C);

	//Generation of latent class labels
	_add_parameter("count_zmisCMK",  CPar_Data_Type::T_INT, &count_zmisCMK, 3, _dat->C, _dat->ncells_mis, K);
	_add_parameter("count_zstratCMK",  CPar_Data_Type::T_INT, &count_zstratCMK, 3, _dat->C, _dat->ncells_strat, K);
	_add_parameter("count_z0CK", CPar_Data_Type::T_INT, &count_z0CK, 2, _dat->C, K);

	//Tabulation
	_add_parameter("count_CK", CPar_Data_Type::T_INT, &count_CK, 2, _dat->C, K);
	_add_parameter("count_CJK2", CPar_Data_Type::T_INT, &count_CJK2, 4, _dat->C, _dat->J - 1, K, 2); 

	//Generation of stratified unobserved counts (0-cells)
	_add_parameter("n0C", CPar_Data_Type::T_INT, &n0C, 1, _dat->C);
	_add_parameter("prob_zeroC", CPar_Data_Type::T_DOUBLE, &prob_zeroC, 1, _dat->C);

	//Parameters
	_add_parameter("log_lambdaCJK2", CPar_Data_Type::T_DOUBLE, &log_lambdaCJK2, 4, _dat->C, _dat->J - 1, K, 2);
	_add_parameter("log_nuCK", CPar_Data_Type::T_DOUBLE, &log_nuCK, 2, _dat->C, K);
	_add_parameter("alphaC", CPar_Data_Type::T_DOUBLE, &alphaC, 1, _dat->C);
	_add_parameter("rhoC", CPar_Data_Type::T_DOUBLE, &rhoC, 1, _dat->C);

	//hyperparameters
	_add_existing_scalar("K", CPar_Data_Type::T_INT, &K);
	_add_existing_scalar("N_Max_factor", CPar_Data_Type::T_DOUBLE, &N_Max_factor);
	_add_existing_scalar("min_n", CPar_Data_Type::T_INT, &min_n);
	_add_existing_scalar("min_lists", CPar_Data_Type::T_INT, &min_lists);
	_add_existing_scalar("a_alpha", CPar_Data_Type::T_DOUBLE, &a_alpha);
	_add_existing_scalar("b_alpha", CPar_Data_Type::T_DOUBLE, &b_alpha);
	_add_existing_scalar("a_lambda", CPar_Data_Type::T_DOUBLE, &a_lambda);
	_add_existing_scalar("b_lambda", CPar_Data_Type::T_DOUBLE, &b_lambda);

	//
	determine_active_lists(_dat, min_n, min_lists);
}

void CParams_NPLCM_CR_Strat_Freq::determine_active_lists(
	CData_DM_Strat* _dat, const int& min_n, const int& min_lists
){
	//Determine the lists that have data for each stratum
	(*this)["activeCJ"].fill(0);
	for (int m = 0; m < _dat->ncells_strat; m++){
		int* X_m = _dat->cells_strat[m] + 1;
		int F_m = _dat->freqM[m];
		int c = _dat->cells_strat[m][0];
		for (int j = 0; j < _dat->J - 1; j++){
			activeCJ[c][j] += X_m[j] == 1 ? F_m : 0;
		}
	}
	for (int c = 0; c < _dat->C; c++){
		for (int j = 0; j < _dat->J - 1; j++){
			activeCJ[c][j] = activeCJ[c][j] >= min_n;
		}
		//for debugging
		std::vector<int> aJ(activeCJ[c], activeCJ[c] + _dat->J - 1);
		//aJ;
		int nl = std::accumulate(aJ.begin(), aJ.end(), 0);
		activeC[c] = nl >= min_lists;
	}

}


void CParams_NPLCM_CR_Strat_Freq::initizalize(gsl_rng *r){
	const int C = (*this)["n0C"].get_dim_lengths()[0];

	//intialize parameters at reasonable values.
	(*this)["count_zmisCMK"].fill(0); 
	(*this)["count_zstratCMK"].fill(0); 
	(*this)["count_z0CK"].fill(0);

	(*this)["count_CK"].fill(0);
	(*this)["count_CJK2"].fill(0);

	(*this)["n_misMC"].fill(0);
	(*this)["n_misC"].fill(0);

	//(*this)["prob_zeroC"].fill(0.0);
	(*this)["n0C"].fill(0);

	(*this)["log_lambdaCJK2"].fill(log(1.0 - 1e-2));
	(*this)["alphaC"].fill(1); //initialize at full spread: uniform prior.
	//Verify this code
	(*this)["rhoC"].fill(1.0 / double(C)); 
}


//--------------------------------------------------------------------------------
// IMPLEMENTATION OF CLASS CNPLCM_CR_Strat_Freq
//--------------------------------------------------------------------------------

//--------------------------------------------------------
// 0) EXECUTION OF GIBBS ITERATION

void CNPLCM_CR_Strat_Freq::Update() {
	reset_counters();

	//latent classification for observed individuals
	sam_zstratCMK();

	//impute missing strata
	sam_n_misMC(); 
	sam_zmisCMK();

	//	imputation of unbserved individuals
	compute_probs_miss();
	try{
		sam_n0C();
	} catch( std::exception& e){
		throw;
	}
	sam_z0CK();

	// Structural parameters
	sam_rhoC();
	sam_lambdaCJK();
	sam_nuCK();
	sam_alphaC();
}

//--------------------------------------------------------
// 1) STRUCTURAL PARAMETERS ------------------------------
//--------------------------------------------------------
inline
void CNPLCM_CR_Strat_Freq::sam_lambdaCJK(){
	const int ncols = data->J - 1;
	(*par)["log_lambdaCJK2"].fill(NAN);
	for(int c = 0; c < data->C; c++){
		if (!par->activeC[c]) continue;
		for(int j = 0; j < ncols; j++){
			if(!par->activeCJ[c][j]) continue;
			for(int k = 0; k < par->K; k++){
				//add prior and cast and sample!
				double a =	par->a_lambda + double(par->count_CJK2[c][j][k][1]);
				double b =	par->b_lambda + double(par->count_CJK2[c][j][k][0]);
				double A = dan_log_gamma_1(r, a);
				double B = dan_log_gamma_1(r, b);
				double C = dan_log_sum(A, B);
				par->log_lambdaCJK2[c][j][k][1] = A - C;
				par->log_lambdaCJK2[c][j][k][0] = B - C;

			}
		}
	}
}

inline
void CNPLCM_CR_Strat_Freq::sam_nuCK(){
	double b = 0.0, a = 0.0, lgamma1, lgamma2, lsumgamma;
	for(int c = 0; c < data->C; c++){
		if (!par->activeC[c]) continue;
		int n_acc = 0;
		double l_acc_prod = 0.0;
		int n_c = data->nC[c] + par->n_misC[c] + par->n0C[c];
		for (int k = 0; k < par->K - 1; ++k){
			int n_ck = par->count_CK[c][k];
			n_acc += n_ck;
			a = double(1 + n_ck);
			b = par->alphaC[c] + double(n_c - n_acc);
			lgamma1 = dan_log_gamma_1(r, a);
			lgamma2 = dan_log_gamma_1(r, b);
			lsumgamma = dan_log_sum(lgamma1, lgamma2);
			par->log_nuCK[c][k] = (lgamma1 - lsumgamma + l_acc_prod);
			l_acc_prod += lgamma2 - lsumgamma; //log(1-V)
		}
		par->log_nuCK[c][par->K -1] = l_acc_prod;
	}
}


inline
void CNPLCM_CR_Strat_Freq::sam_alphaC(){
	//update alphaC[]
	for (int c = 0; c < data->C; c++){
		if (!par->activeC[c]) continue;
		par->alphaC[c] = gsl_ran_gamma(
			r, par->a_alpha + double(par->K) - 1.0, 
			1.0 / (par->b_alpha - par->log_nuCK[c][par->K - 1])
		);
	}
}

//--------------------------------------------------------
// 2) LATENT CLASS LABELS   ------------------------------
//        
inline 
void CNPLCM_CR_Strat_Freq::reset_counters(){
	//Initialize all tabulation holders:
	(*par)["count_CK"].fill(0);
	(*par)["count_CJK2"].fill(0);
}

inline
void CNPLCM_CR_Strat_Freq::sam_z0CK(){
	//needs: n0C, lambdas, nus, 
	const int ncols = data->J -1;
	CVariable_Container& l_lambdaCJK2 = (*par)["log_lambdaCJK2"];
	int zero = 0;

	for(int c = 0; c < data->C; c++){
		std::vector<double> table_z(par->K, 0);
		if (!par->activeC[c]) continue;
		if (par->n0C[c] < 1) {
			std::fill(par->count_z0CK[c], par->count_z0CK[c] + par->K, 0);
			continue;
		}
		int* dest_c =  par->count_z0CK[c];
		double lg = -INFINITY;
		for (int k = 0; k < par->K; k++){
			double lp = par->log_nuCK[c][k];
			for(int j = 0; j < ncols; j++){
				if (!par->activeCJ[c][j]) continue;
				//lp += par->log_lambdaCJK2[c][j][k][0];
				lp += l_lambdaCJK2.arr_elem<double>(c, j, k, zero);
			}
			table_z[k] = lp;
			lg = lg < lp ? lp : lg;
		}
		for (int k = 0; k < par->K; k++) {
			table_z[k] = exp(table_z[k] - lg);
		}
		gsl_ran_multinomial(r, par->K, par->n0C[c], &(table_z[0]), (unsigned int*)(dest_c));

		//update counters
		for (int k = 0; k < par->K; k++){
			par->count_CK[c][k] += dest_c[k];
			for (int j = 0; j < ncols; j++){
				if(!par->activeCJ[c][j]) continue;
				par->count_CJK2[c][j][k][0] += dest_c[k];
			}
		}
	}
}

inline
void CNPLCM_CR_Strat_Freq::sam_zstratCMK(){
	//update count_zstratCMK[][][], count_CK[][], count_CJK2[][][][]
	std::vector<double> probs(par->K, 0);
	const int ncols = data->J - 1;
	CVariable_Container& l_lambdaCJK2 = (*par)["log_lambdaCJK2"];

	(*par)["count_zstratCMK"].fill(0);
	for (int m = 0; m < data->ncells_strat; m++){
		int  c = data->cells_strat[m][0];
		if (!par->activeC[c]) continue;
		int* x_m = data->cells_strat[m] + 1;
		int  F_m = data->freqM_strat[m];
		int* dest_cm = par->count_zstratCMK[c][m];
		double lg = -INFINITY;
		if (!_is_zeros_strat(x_m, c)){
			for (int k = 0; k < par->K; k++){
				double lprod = par->log_nuCK[c][k];
				for (int j = 0; j < ncols; j++){
					if (!par->activeCJ[c][j]) continue;
					//lprod += par->log_lambdaCJK2[c][j][k][ x_m[j] ];
					lprod += l_lambdaCJK2.arr_elem<double>(c, j, k, x_m[j]);
				}
				lg = lg < lprod ? lprod : lg;
				probs[k] = lprod;
			}
			for (int k = 0; k < par->K; k++){
				probs[k] = exp(probs[k] - lg);
			}
			//sample!
			gsl_ran_multinomial(r, par->K, F_m, &(probs[0]), (unsigned int*) dest_cm);
			//update counters
			for (int k = 0; k < par->K; k++){
				par->count_CK[c][k] += dest_cm[k];
				for (int j = 0; j < ncols; j++){
					if (!par->activeCJ[c][j]) continue;
					par->count_CJK2[c][j][k][ x_m[j] ] += dest_cm[k];
				}
			}

		} else {
			//count the extra n0c
			// COMPLETE THIS CODE
			//std::vector<int> Z(x_m, x_m + ncols);
			//int n0c = F_m;
		}
	}
}

inline
void CNPLCM_CR_Strat_Freq::sam_zmisCMK(){
	//update count_zmisCMK[][][], count_CK[][], count_CJK2[][][][]
	std::vector<double> probs(par->K, 0);
	const int ncols = data->J - 1;
	(*par)["count_zmisCMK"].fill(0);
	CVariable_Container& l_lambdaCJK2 = (*par)["log_lambdaCJK2"];


	for (int m = 0; m < data->ncells_mis; m++){
		int* x_m = data->cells_mis[m] + 1;
		for (int c = 0; c < data->C; c++){
			if (!par->activeC[c]) continue;
			int  F_mc = par->n_misMC[m][c];
			int* dest_mc = par->count_zmisCMK[c][m];
			if (F_mc == 0) {
				std::fill(dest_mc, dest_mc + par->K, 0);
			} else if (_is_zeros_strat(x_m, c)) {
				//add to extra n0s.
				// COMPLETE TIS CODE
				//std::vector<int> v(x_m, x_m + ncols);
				//F_mc;
			} else {
				double lg = -INFINITY;
				for (int k = 0; k < par->K; k++){
					double lprod = par->log_nuCK[c][k];
					for (int j = 0; j < ncols; j++){
						if (!par->activeCJ[c][j]) continue;
						//lprod += par->log_lambdaCJK2[c][j][k][ x_m[j] ];
						lprod += l_lambdaCJK2.arr_elem<double>(c, j, k, x_m[j]);
					}
					lg = lg < lprod ? lprod : lg;
					probs[k] = lprod;
				}
				for (int k = 0; k < par->K; k++){
					probs[k] = exp(probs[k] - lg);
				}
				//sample!
				gsl_ran_multinomial(r, par->K, F_mc, &(probs[0]), (unsigned int*) dest_mc);
				//update counters
				for (int k = 0; k < par->K; k++){
					par->count_CK[c][k] += dest_mc[k];
					for (int j = 0; j < ncols; j++){
						if(!par->activeCJ[c][j]) continue;
						par->count_CJK2[c][j][k][ x_m[j] ] += dest_mc[k];
					}
				}
			}
		}
	}
}

//--------------------------------------------------------
// 3) UNOBSERVED POPULATION SAMPLING ---------------------
inline
void CNPLCM_CR_Strat_Freq::compute_probs_miss(){
	//update prob_zeroC[]
	const int ncols = data->J - 1;
	std::vector<double> terms_cK(par->K, 0.0);

	for (int c = 0; c < data->C; c++){
		par->prob_zeroC[c] = 0.0;
		if (!par->activeC[c]) continue;
		double lmax_prob = -INFINITY;
		int index_max = -1;
		for (int k = 0; k < par->K; ++k) {
			double lprod = par->log_nuCK[c][k];
			for (int j = 0; j < ncols; j++) {
				if (!par->activeCJ[c][j]) continue;
				lprod += par->log_lambdaCJK2[c][j][k][0];
			}
			if (lmax_prob < lprod){
				lmax_prob = lprod;
				index_max = k;
			}
			terms_cK[k] = lprod;
		}
		double sum = 0.0;
		for (int k = 0; k < par->K; k++){
			if (k == index_max) continue;
			sum += exp(terms_cK[k] - lmax_prob);
		}
		double logprob = lmax_prob + log1p(sum);
		par->prob_zeroC[c] = par->rhoC[c] * exp(logprob);
	}
}

	//
	//
	inline void CNPLCM_CR_Strat_Freq::sam_n0C()
	{
		//Sample from negative multinomial.
		double p0 = 0.0;
		int n = 0;
		unsigned int n0 = 0;

		for (int c = 0; c < data->C; c++){
			if (!par->activeC[c]) continue;
			p0 += par->prob_zeroC[c];
			n += data->nC[c] + par->n_misC[c];
		}
		int counter = 0;
		do{
			if (counter++ > 100000)
				throw(std::runtime_error("Error: p0 too large"));
			n0 = gsl_ran_negative_binomial(r, 1.0 - p0, n);
		} while (n0 > (par->N_Max_factor - 1) * n);
		gsl_ran_multinomial(r, data->C, n0, par->prob_zeroC, (unsigned int *)par->n0C);
		//dan_neg_multinomial(r, n, data->C, par->prob_zeroC, par->n0C);
}


//--------------------------------------------------------
// 4) IMPUTATION OF MISSING STRATIFICATION ---------------

inline
void CNPLCM_CR_Strat_Freq::sam_n_misMC(){
	//Imputation of missing strata
	//Here we need to set par->n_misMC and par->n_misC
	(*par)["n_misMC"].fill(0);
	(*par)["n_misC"].fill(0);
	CVariable_Container& l_lambdaCJK2 = (*par)["log_lambdaCJK2"];

	for (int m = 0; m < data->ncells_mis; m++){
		std::vector<double> prob_mC(data->C, 0.0);
		int* X_m = data->cells_mis[m] + 1;// cells_mis[m][0] is -1
		int n_m = data->freqM_mis[m];
		int* dest_mC = par->n_misMC[m];
		std::fill(prob_mC.begin(), prob_mC.end(), 0.0);
		for (int c = 0; c < data->C; c++){
			if (!par->activeC[c]) continue;
			if (_is_zeros_strat(X_m, c)){
				prob_mC[c] = 0.0;
			} else {
				double sum = 0.0;
				for (int k = 0; k < par->K; k++){
					double lp1 = par->log_nuCK[c][k];
					for (int j = 0; j < data->J - 1; j++){ 
						if (!par->activeCJ[c][j]) continue;
						//lp1 += par->log_lambdaCJK2[c][j][k][ X_m[j] ];
						lp1 += l_lambdaCJK2.arr_elem<double>(c, j, k, X_m[j]);
					}
					sum += exp(lp1);
				}
				prob_mC[c] = par->rhoC[c] * sum;
			}
		}
		gsl_ran_multinomial(r, data->C, n_m, &(prob_mC[0]), (unsigned int *)dest_mC);
		//count
		for (int c = 0; c < data->C; c++){
			par->n_misC[c] += dest_mC[c];
		}
	}
}

inline
void CNPLCM_CR_Strat_Freq::sam_rhoC(){
	// updates rhoC[]
	std::vector<double> aC(data->C, 1.0); //1.0 is flat prior distribution

	for (int c = 0; c < data->C; c++){
		if (!par->activeC[c]) {
			aC[c] = 0.0;
			continue;
		}
		for (int k = 0; k < par->K; k++){
			aC[c] += par->count_CK[c][k];
		}
	}
	gsl_ran_dirichlet(r, data->C, &(aC[0]), par->rhoC);
}


void CNPLCM_CR_Strat_Freq::Initializes(){
	//Let's use LCMCRs on observed data for initialization
	const int ncols = data->J - 1;
	std::vector<std::vector<int> > XC(data->C);
	std::vector<int> n_activC(data->C);

	//read data and put it into data->C flat arrays (RMO)
	for (int c = 0; c < data->C; c++){
		n_activC[c] = std::accumulate(par->activeCJ[c], par->activeCJ[c] + ncols, 0);
		XC[c].resize(data->nC[c] * n_activC[c]);
	}

	std::vector<int> indexC(data->C, 0);
	for(int m = 0; m < data->ncells_strat; m++){
		int c = data->cells_strat[m][0];
		int* X_mJ = data->cells_strat[m] + 1;
		int F_m = data->freqM_strat[m];
		if(!par->activeC[c]) continue;
		if(_is_zeros_strat(data->cells_strat[m] + 1, c)){
			continue;
		}
		for (int i = 0; i < F_m; i++){
			
			for (int j = 0, jj = 0; j < ncols; j++){
				if (par->activeCJ[c][j]){
					XC[ c ][ indexC[c]*n_activC[c] + jj++] = X_mJ[j];
				}
			}
			indexC[c]++;
		}
	}

	//Create C basic models and initialize them.
	std::vector<CData_DM*> data_s(data->C);
	std::vector<CNPLCM_CR_Basic_Freq*> model_s(data->C);
	for (int c = 0; c < data->C; c++){
		if(!par->activeC[c]) continue;
		CData_DM* d = new CData_DM();
		d->Set_Manually(&(XC[c][0]), n_activC[c], indexC[c], data->levelsJ + 1);
		CNPLCM_CR_Basic_Freq* m = new CNPLCM_CR_Basic_Freq(d, par->K, par->a_alpha, par->b_alpha);
		m->reseed_rng(this->ran_seed);
		m->Initializes();
		model_s[c] = m;
		data_s[c] = d;
	}

	//Copy the initialized values to strata. Don't need to copy any counts with K.
	(*par)["n_misC"].fill(0);
	(*par)["n_misMC"].fill(0);
	int act;
	act = std::accumulate(par->activeC, par->activeC + data->C, 0);
	for (int m = 0; m < data->ncells_mis; m++){
		//just split in K equal parts
		int nc  = data->freqM_mis[m] / act; //integer division
		int rem = data->freqM_mis[m] % act;
		int cc = 0;
		for (int c = 0; c < data->C; c++){
			if(!par->activeC[c]) continue;
			par->n_misC[c] += par->n_misMC[m][c] = nc + (cc++ == 0 ? rem : 0);
		}
	}

	(*par)["log_lambdaCJK2"].fill(NAN);
	(*par)["log_nuCK"].fill(NAN);
	for (int c = 0; c < data->C; c++){
		if (!par->activeC[c]){
			par->alphaC[c] 		= 1.0; 
			par->prob_zeroC[c]  = 0.0;
			par->n0C[c] 		= 0.0; 
			par->rhoC[c] 		= 0.0; 
		}  else {
			CParams_NPLCM_CR_Basic_Freq* p = model_s[c]->par;
			//CData_DM* d = data_s[c];
			par->alphaC[c] 		= p->alpha;
			par->prob_zeroC[c]  = p->prob_zero;
			par->n0C[c] 		= p->n0; 
			int jj = 0;
			for(int j = 0; j < ncols; j++){
				if (!par->activeCJ[c][j]) continue;
				for(int k = 0; k < par->K; k++){
					par->log_lambdaCJK2[c][j][k][0] = p->log_lambdaJK2[jj][k][0];
					par->log_lambdaCJK2[c][j][k][1] = p->log_lambdaJK2[jj][k][1];
				}
				jj++;
			}
			std::copy(p->log_nuK, p->log_nuK + par->K, par->log_nuCK[c]);
			std::vector<double> v(p->log_nuK, p->log_nuK + par->K);
			par->rhoC[c] = 1.0 / double(act);
		}
	}

	//Deallocate constructed objects
	for (int c = 0; c < data->C; c++){
		if (!par->activeC[c]) continue;
		delete model_s[c];
		delete data_s[c];
	}
}
