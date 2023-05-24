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

#ifndef _DAN_MATH_GSL_H
#define _DAN_MATH_GSL_H

//Mscelaneous math functions that depend on GSL.

#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include "dan_sys.h"
#include "dan_math.h"


//Prototypes
gsl_rng *dan_initRandom();
double dan_log_sum(const double &a, const double &b);
double dan_log_sum_K(double * restrict x, const int& K);
double dan_inv_digamma(const double& y);
template<bool LOG_SCALE> double dan_normalize(double *src, double *dst, const int& K);
double dan_beta1alpha_trn(gsl_rng* r, double alpha, double a,double b, bool log_one_minus_beta = false);
double dan_log_gamma_1(const gsl_rng* r, const double &shape);
void dan_log_dirichlet(const gsl_rng* r, const int &K, const double * restrict alpha, double * restrict theta);
double dan_beta1alpha_trn(gsl_rng* r, double alpha, double a,double b, bool log_one_minus_beta);
int dan_multinomial_1(const gsl_rng* r, const int &K, const double * restrict p, bool is_norm);
int dan_multinomial_1_norm(const gsl_rng* r, const int &K, const double * restrict p, const double &norm);
int dan_multinomial_1_cum_prob(const gsl_rng* r, const int &K, const double * restrict cum_p);
double dan_weib(const double &x, const double &alpha, const double &beta);
double dan_cdf_weib(const double &x, const double &alpha, const double &beta);
double dan_inv_cdf_weib(const double &y, const double &alpha, const double &beta);
double dan_lTruncWeibull(const gsl_rng* r, const double &alpha, const double &beta, const double &l_trn);
double dan_lTruncNormal(const gsl_rng* r, const double &l_trn);
void dan_neg_multinomial(const gsl_rng* r, const int &k, const int &size, const double p[], int ns[]);
double dan_lratio_gamma(const double &a, const double &b);
double dan_log1exp(const double &x);
double dan_log_gaussian_cdf(const double &x);



//special functions
inline
double dan_inv_digamma(const double& y){
	//method from Minka, 2000
	double x, xold;
	x = y > -2.22? exp(y) + 0.5 : - 1.0 / (y + 0.577215664901532); //the constant is -psi(1)
	do{
		xold = x;
		x -= (gsl_sf_psi(x) - y) / gsl_sf_psi_1(x);
	}while(fabs((x-xold)/x)> 1e-10);
	return(x);
}

inline 
double dan_log_gaussian_cdf(const double &x){
	//approximation for small values between -1000 and -30.
	return (x > -30.0 ?  log(gsl_cdf_ugaussian_P(x))
		: -4.20089235997904 +
			x * (  0.0154449653185494 +
			x * (-0.499966366781208 +
			x * ( 3.59656615605398e-08 + 
			x * 1.4219631098393e-11 )))
			);

}


//Numerical tricks
inline 
double dan_log_sum(const double &a, const double &b){
	// = log( exp(a) + exp(b))
	return a < b ? b + gsl_sf_log_1plusx(exp(a - b)) : a + gsl_sf_log_1plusx(exp(b - a));
}

inline
double dan_log_sum_K(double * restrict x, const int& K){
	// = \log( \sum_i \exp(x_i))
	double res = -DAN_INF;
	double * restrict ptr = x;
	for (int k = 0; k < K; k++) res = dan_log_sum(res, *(ptr++));
	return(res);
}

//analyze this. Consider cases: x<log(DLB_EPSILON), log(DBL_EPSILON) < x < -log(DBL_+EPSILON), x>max_exp_arg
//inline
//double dan_log1exp(const double &x){
//	return (x < MAX_EXP) ? 
//		log( 1.0 + exp(x) ) : x + gsl_sf_log_1plusx(exp(-x));
//}
inline
double dan_log1exp(const double &x){
	return (x < MAX_EXP) ? 	
			gsl_sf_log_1plusx(x) : (
				(x < -MAX_EXP) ? 
				log(1.0 + exp(x)) : 
				x + gsl_sf_log_1plusx(exp(-x))
			);
}


inline
double dan_lratio_gamma(const double &a, const double &b){
 return ((a<1e-300)|(b<1e-300)) ?
		gsl_sf_lngamma(a + 1.0) - gsl_sf_lngamma(b + 1.0) - log(a) + log(b):
		gsl_sf_lngamma(a)-gsl_sf_lngamma(b);

}
template<bool LOG_SCALE>
inline double dan_normalize(double *src, double *dst, const int& K){
	double norm = 0.0;
	int k;
	if (LOG_SCALE){
		norm  = dan_log_sum_K(src, K); //initial value of norm doesn't matter here.
		for (k = 0; k < K; k++) dst[k] = src[k] - norm;
	} else {
		for (k = 1; k < K; k++) norm += src[k];
		for (k = 0; k < K; k++) dst[k] = src[k]/norm;
	}
	return(norm);
}
///////////////////////////////////////////////////////////
/* 
	Numeros aleatorios.
	Rutinas para hacer las cosas un poco mas faciles con las 
	 funciones de generacion de numeros aleatorios de gsl.
*/

//Inicializa un generador de numeros aleatorios.

inline
gsl_rng* dan_initRandom(){
	/* create and initializes a rng object and
	seeds the random number generator with the system time
	*/
	gsl_rng *r;
	r = gsl_rng_alloc(gsl_rng_taus2);
	int off = sizeof(time_t) / 2;
	int t = ( int(time(0)) << off) >> off ;
	gsl_rng_set(r,t);
	return(r);
}

inline
gsl_rng* dan_initRandom(int seed){
	/* create and initializes a rng object and
	seeds the random number generator with the system time
	*/
	gsl_rng *r;
	r = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(r,seed);
	return(r);
}

inline
void dan_log_dirichlet(const gsl_rng* r, const int &K, const double * restrict alpha, double * restrict theta){
	int  i;
	double norm = 0.0;

	for (i = 0; i < K; i++) theta[i] = dan_log_gamma_1(r, alpha[i]);
	double max_alpha = dan_max_vector_double(theta, K);
	for (i = 0; i < K; i++) norm += exp(theta[i] - max_alpha);
	norm = log(norm) + max_alpha;
	for (i = 0; i < K; i++) theta[i] = theta[i] - norm;
}

inline
double dan_log_gamma_1(const gsl_rng* r, const double &shape){
	//function to sample from y such that e^y ~ gamma (shape = alpha, scale = 1)
	// for small values of shape uses the method from p.45 of RObert & Casella (2002)
	return (
		(shape < 0.5) ?
			log(gsl_rng_uniform_pos (r))/shape + log(gsl_ran_gamma(r, 1.0 + shape, 1.0))
			:	log(gsl_ran_gamma(r, shape, 1.0))
	);// result;
}


inline
int dan_multinomial_1(const gsl_rng* r, const int &K, const double * restrict p, bool is_norm){
	double norm = 0.0, cum[10000];
	int k;
	for (k = 0; k < K; ++k) cum[k] = norm += p[k]; 
	double u = gsl_rng_uniform(r) * norm; //a number between 0 and norm.
	for(k = 0; cum[k] <= u; ++k);
	return k;
}
inline
int dan_multinomial_1_cum_prob(const gsl_rng* r, const int &K, const double * restrict cum_p){
	int k;
	double u = gsl_rng_uniform(r) * cum_p[K-1]; //a number between 0 and norm.
	k = 0;
	while(u > cum_p[k]) k++;
	return k;
}
inline
int dan_multinomial_1_norm(const gsl_rng* r, const int &K, const double * restrict p, const double &norm){
	int k=0;
	double u = gsl_rng_uniform(r) * norm; //a number between 0 and norm.
	//A simple linear search. this should be good enough for low values of K.
	double cum = p[0];
	while(u > cum){
		cum += p[++k];
	}
	return k;
}
inline 
double dan_beta1alpha_trn(gsl_rng* r, double alpha, double a,double b, bool log_one_minus_beta){
	//samples from truncated (to a<b) beta(1,alpha) distribution.
	// If log_one_minus_beta = true, it returns log(1-beta), with beta ~ Beta(beta|1,alpha)*1(a<beta<b)
	double u = gsl_ran_flat(r,0.0, 1.0);
	//New version
	double a1 = log(u) + alpha * log(1.0 - b);
	double a2 = log(1.0 - u) + alpha * log(1.0 - a);
	double tmp = dan_log_sum(a1, a2) / alpha;
	double res = log_one_minus_beta ? tmp : 1.0 - exp(tmp);
	return(res);
}
inline
void dan_neg_multinomial(const gsl_rng* r, const int &k, const int &size, const double p[], int ns[]){
	//k: count in conditioning cell
	//size: number of random cells
	//p: probabilities of random cells (note sum p < 1)
	//ns: k returning values.
	//for k = 1 this should be a regular negative binomial.
	double norm = 0.0;

	for (int i=0; i < size; i++){
		norm += p[i];
	}
	double gamma = gsl_ran_gamma(r, (double) k, 1.0);
	for (int i = 0; i < size; i ++){
		ns[i] = gsl_ran_poisson(r, gamma * p[i]/(1.0-norm));
	}
	return;
}
inline
int dan_neg_binomial_p0(const gsl_rng* r, const int &size, const double p0){
	// Number of failures before (size) successes
	// Sece
	//p: probability of success

	double gamma = gsl_ran_gamma(r, double(size), 1.0);
	double f = exp(log(p0)- log1p(-p0));
	return gsl_ran_poisson(r, gamma * f);
}
inline
double dan_lTruncWeibull(const gsl_rng* r, const double &alpha, const double &beta, const double &l_trn){
	double u = gsl_ran_flat(r, dan_cdf_weib(beta*l_trn, alpha, 1.0), 1.0);//beta = 1.0. Then scale.
	return(dan_inv_cdf_weib(u, alpha, 1.0)/beta);
}

inline
double dan_cdf_weib(const double &x, const double &alpha, const double &beta){
	return( 1.0-exp(-pow(beta*x,alpha)));
}
inline
double dan_weib(const double &x, const double &alpha, const double &beta){
	return (	pow(alpha*beta,alpha)*pow(x,alpha-1.0)*exp(-pow(x*beta,alpha)));
}
inline
double dan_inv_cdf_weib(const double &y, const double &alpha, const double &beta){
	return (pow(-log(1.0-y),1.0/alpha)/beta);
}

inline
double dan_lTruncNormal(const gsl_rng* r, const double &l_trn){
	//Adapted from method from Robert 1995
	double z, alpha;
	if (l_trn <= 0.0){
		do{
			z = gsl_ran_ugaussian(r);
		} while (z < l_trn);	
	} else {
		alpha = 0.5 * (l_trn + sqrt(l_trn * l_trn + 4.0));
		do{
			z = -log(gsl_rng_uniform_pos(r))/alpha + l_trn; //sample from translated ~exp(alpha)
		} while ( log(gsl_rng_uniform_pos(r)) > -0.5 * (z - alpha)*(z - alpha));
	} 
	return(z);
}
#endif
