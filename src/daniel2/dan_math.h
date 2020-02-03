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

#ifndef _DAN_MATH_H
#define _DAN_MATH_H

//Miscelaneous math functions that don't depend on GSL.
#include <math.h>
#include <float.h>
#include "dan_sys.h"

//Numerical precision constants
#define MAX_EXP log(DBL_EPSILON) //max x such that exp(x) + 1 gets rounded to 1
#define MAX_EXP_ARG log(DBL_MAX) // Max x such that exp(x) < inf
#define MIN_EXP_ARG log(-DBL_MAX) // Min x such that exp(x) > 0
#ifndef DAN_INF
#define DAN_INF HUGE_VAL
#endif

//prototypes.
double logit_fn(const double &x);
double alogit_fn(const double &x);
template<typename T> T dan_max(T *vec, const int &length);
template<typename T> T dan_min(T *vec, const int &length);
template<typename T> int dan_argmax(T *vec, const int &length);
template<typename T> int dan_argmin(T *vec, const int &length);
template<typename T> T dan_prod_vec(T *data, int &n);
template<typename T> T dan_sum_vec(T *data, int &n);

//A simple template iterator. Applies "func" to each element of src and places the results on dst.
// example: dan_apply<double,double,log>(x, y, 10);
template<typename T_src, typename T_dst, T_dst(&func)(T_src)> void dan_apply(T_src* src, T_dst* dst, const int& length);


inline
double alogit_fn(const double &x){
	//return exp(x)/(1.0+exp(x));
	return x > 1 ? 1.0 / (1.0 + exp(-x)) : 1.0 - 1.0 / (1.0 + exp(x));
}
inline
double logit_fn(const double &x){
	return -log(1.0 / (1.0/x - 1.0));
}

template<typename T>
inline int dan_argmin(T *vec, const int &length){
	int argmin = 0;
	T min = vec[0];
	for (int i = 1; i < length; i++){
		if (vec[i] < min) {
			min = vec[i];
			argmin = i;
		}
	}
	return(argmin);
}

template<typename T>
inline int dan_argmax(T *vec, const int &length){
	int argmax = 0;
	T max;
	max = vec[0]; 
	for (int i = 1; i < length; i++){
		if (vec[i] > max) {
			max = vec[i];
			argmax = i;
		}
	}
	return(argmax);
}
template<typename T>
inline T dan_min(T *vec, const int &length){
	T min = vec[0];
	for (int i = 1; i < length; i++){
		if (vec[i] < min) min = vec[i];
	}
	return(min);
}
template<typename T>
inline T dan_max(T *vec, const int &length){
	T max;
	max = vec[0]; 
	for (int i = 1; i < length; i++){
		if (vec[i] > max) max = vec[i];
	}
	return(max);
}
template<typename T>
inline T dan_sum_vec(T * restrict data, int &n){
	T p = (T)0;
	for (int i = 0; i < n; i++){
		p += data[i];
	}
	return(p);
}
template<typename T>
inline T dan_prod_vec(T * restrict data, int &n){
	T p = (T)1;
	for (int i = 0; i < n; i++){
		p *= data[i];
	}
	return(p);
}

template<typename T_src, typename T_dst, T_dst(&func)(T_src)>
inline void dan_apply(T_src* src, T_dst* dst, const int& length){
	for (int i = 0; i < length; i ++){
		dst[i] = func(src[i]);
	}
}
inline
double dan_Kahns_sum(double *x, const int& length){
	//From wikipedia:
	//function KahanSum(input)
 //   var sum = 0.0
 //   var c = 0.0          //A running compensation for lost low-order bits.
 //   for i = 1 to input.length do
 //       var y = input[i] - c    //So far, so good: c is zero.
 //       var t = sum + y         //Alas, sum is big, y small, so low-order digits of y are lost.
 //       c = (t - sum) - y   //(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
 //       sum = t             //Algebraically, c should always be zero. Beware overly-aggressive optimising compilers!
 //       //Next time around, the lost low part will be added to y in a fresh attempt.
 //   return sum
	double sum = 0.0, c = 0.0, t, y;
	for (int i = 0; i < length; i ++){
		y = x[i] - c;
		t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	return (sum);
}

#define dan_log_sum_exps dan_log_sum_exps2
inline
double dan_log_sum_exps1(double* x, const int& length){
	//This uses Kahn's summation:
	double res = 0.0, c = 0.0, y, t;
	double m = dan_max<double>(x, length); //normalizing factor for EXP
	for (int i = 0; i < length; i++){
		y = exp(x[i] - m) - c;
		t = res + y;
		c = (t - res) - y;
		res = t;
	}
	return(log(res) + m);
}

inline
double dan_log_sum_exps2(double* x, const int& length){
	//This uses Kahn's summation:
	double res = 0.0, c = 0.0, y, t;
	double m = dan_max<double>(x, length); //normalizing factor for EXP
	for (int i = 0; i < length; i++){
		y = exp(x[i] - m) - c;
		t = res + y;
		c = (t - res) - y;
		res = t;
	}
	return(log(res) + m);
}

//////////////////////////////////////////////////////////////////////////
//DEPRECATED STUFF. 


//Use the templates instead.

inline int dan_argmin_vector_double(double *vec, const int &length){return dan_argmin<double>(vec, length);}
inline int dan_argmax_vector_int(int *vec, const int &length){return dan_argmax<int>(vec, length);}
inline int dan_max_vector_int(int *vec, const int &length){return dan_max<int>(vec, length);}
inline double dan_min_vector_double(double const vec[], const int &length){return dan_min<double>(const_cast<double*>(vec), length);}
inline double dan_max_vector_double(double const vec[], const int &length){return dan_max<double>(const_cast<double*>(vec), length);}


#endif
