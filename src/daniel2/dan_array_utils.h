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

#ifndef _ARRAY_UTILS_H
#define	_ARRAY_UTILS_H

#include <string>
#include <cstring>
#include <stdlib.h>
#include <cstdarg>
#include <vector>

//////////////////////////////////
//VECTOR AND MATRIX UTILITIES.
//  **Remember -> C: rightmost varying fastest R: leftmost varying fastest
//	Takes a flat vector and reshapes it as an array. Rightmost varying the fastest
//	Example:
//		double ****Array = (double****)dan_flat2array(MyVector, sizeof(double), 4, 3,4,5,6);
//		takes MyVector (3x4x5x6 elements) and shapes it into Array[3][4][5][6].
//		Notes:  Both Matrix and MyVector point to the same memory region.
//	The array pointer overhead MUST be cleaned up manually:
//		free(Array);

//Malloc utilities
void *dan_malloc(int bytes, char const * name_var, char const * name_function);
void *dan_malloc(int bytes);
void dan_free(void* ptr);

// compile time type indirector. Returns the _T plus n levels of indirection e.g.
// dan_n_pointer<int,4>::type p; is equivalent to int**** p;
template <typename _T, int n>
class dan_n_pointer {
public:
	typedef typename dan_n_pointer<_T, n - 1>::type* type;
};
template <typename _T>
class dan_n_pointer<_T, 0> {
public:
    typedef _T type;
};


//Creation of multidimensional arrays
// versions with malloc
void *dan_flat2arrayND_ln(void* data, int size_elem, int dims, int *lengths);
void *dan_flat2arrayND(void* data, int size_elem, int dims, ...);
void dan_free_arrayND_overhead(void* arreglo);
// versions with "operator new"
void *dan_flat2arrayND_cpp(void* data, int size_elem, const std::vector<int>& lengths);
void *dan_flat2arrayND_cpp(void* data, int size_elem, int dims, ...);


//Manipulation of array indices
void dan_next_in_seq(int* z, int* lengths, int J);
int dan_index_2_linear_position_cstyle(int* ind, int &J, int *levelsJ);
int dan_index_2_linear_position_fortranstyle(int* ind, int &J, int *levelsJ);
template<typename T> void dan_transpose(const T* orig_flat, T* dst_flat, const int& dims, const int* lengths);
template<typename T> void dan_aperm(const T* orig_flat, T* dst_flat, const int& dims, const int lengths[], const int order[]);
template<typename T> void dan_array2R(std::string& s, T* src_flat, const int& dims, int* lengths, const int perline = 10, bool miss = false, const T& mis_code = -1);

//Contingency tables.
void dan_full_contingency_table(int **dataIJ, int &N, int &J, int* output, int *levelsJ, bool C_style = true, bool clear_output = true);

//-----Implementation of inlines and templates--------

inline
void dan_next_in_seq(int* z, int* lengths, int J){
	//computes the next combination of levels in the array z.
	int digit = 0;
	bool done = false;
	do{
		if(++z[digit] == lengths[digit]){
			z[digit++] = 0;
		} else done = true;
	} while(!done);
}

inline
int dan_index_2_linear_position_fortranstyle(int* ind, int &J, int *levelsJ){
	//computes the linear position of an indexed element of an array (col major order)
	// *levelsJ contains the number of levels of each dimension.
	//*ind contains the indices to translate
	int pos = ind[0], power = 1;
	for (int j = 1; j < J; j ++) pos += (power *= levelsJ[j-1])* ind[j];
	return(pos);
}

inline
int dan_index_2_linear_position_cstyle(int* ind, int &J, int *levelsJ){
	//computes the linear position of an indexed element of an array (row major order)
	// *levelsJ contains the number of levels of each dimension.
	//*ind contains the indices to translate
	int pos = ind[J-1];
	int power = 1;
	for (int j = 2; j <= J; j ++){
		power *= levelsJ[J-j+1];
		pos += power * ind[J-j];
	}
	return(pos);
}
inline 
void dan_transpose_untyped(const void* orig_flat, void* dst_flat, const int& dims, const int* lengths, const int& size_elem){
	// transform a multidimensional array from row-major order to col-major-order
	int x, y;
	int index, multiplier, size = 1;
	for (int i = 0; i < dims; i++){
		size *= lengths[i];
	}
	if (dims == 1) {
		std::memcpy(dst_flat, orig_flat, size_elem * lengths[0]);
		return;
	}
	for (int j = 0; j < size; j ++){
		y = 0; //destination index
		x = j; //source index
		multiplier = size / lengths[dims - 1];
		for (int i = dims - 1; i >= 0 ; i --){
			index = x % lengths[i];
			y += index * multiplier;
			x /= lengths[i];
			if (x == 0) break;
			multiplier /= lengths[i - 1];
		}
		//dst_flat[y] = orig_flat[j];
		std::memcpy( (void*)((char*)dst_flat + size_elem * y), (void*)((char*)orig_flat + size_elem * j), size_elem);
	}
}

template<typename T>
void dan_transpose(const T* orig_flat, T* dst_flat, const int& dims, const int* lengths){
	// transform a multidimensional array from row-major order to col-major-order
	int x, y;
	int index, multiplier, size = 1;
	for (int i = 0; i < dims; i++){
		size *= lengths[i];
	}
	if (dims == 1) {
		std::memcpy(dst_flat, orig_flat, sizeof(T) * lengths[0]);
		return;
	}
	for (int j = 0; j < size; j ++){
		y = 0; //destination index
		x = j; //source index
		multiplier = size / lengths[dims - 1];
		for (int i = dims - 1; i >= 0 ; i --){
			index = x % lengths[i];
			y += index * multiplier;
			x /= lengths[i];
			if (x == 0) break;
			multiplier /= lengths[i - 1];
		}
		dst_flat[y] = orig_flat[j];
	}
}

template<typename T>
void dan_aperm(const T* orig_flat, T* dst_flat, const int& dims, const int lengths[], const int order[]){
	// transpose arbitrary indexes from an array.
	int x, y;
	int index, weights[100], size = 1;
	for (int i = 0; i < dims ; i++){
		size *= lengths[i];
	}
	int prod = 1;
	for (int i = dims - 1; i >=0; i--){
		weights[i] = prod;
		prod *= lengths[order[i]];
	}
	if (dims == 1) {
		std::memcpy(dst_flat, orig_flat, sizeof(T) * lengths[0]);
		return;
	}
	for (int j = 0; j < size; j ++){
		y = 0; //destination index
		x = j; //source index
		for (int i = dims - 1; i >= 0 ; i --){
			index = x % lengths[i];
			x /= lengths[i];
			y += index * weights[order[i]];
			if (x == 0) break;
		}
		dst_flat[y] = orig_flat[j];
	}
}


//Do not use anymore---------------
inline
void dan_fill_array(double* x, double Val, const int& length){
	for (int i = 0; i < length; i++) {
		x[i] = Val;
	}
}


//////////////////////////////////////////////////////////////////////////
//DEPRECATED FUNCTIONS.
//Implemented in the *.cpp file
//Use arrayND instead.
void* dan_flat2array2D(void* data, int size_elem, int dim_i, int dim_j);
void* dan_flat2array3D(void* data, int size_elem, int dim_i, int dim_j, int dim_k);
void dan_free_array2D_overhead(void* myArray);
void dan_free_array3D_overhead(void* myArray,int dim_j);
#endif
