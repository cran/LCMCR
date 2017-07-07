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

#include "dan_array_utils.h"
#include "dan_sys.h"
/*----------------------------------------------------------
	Arrays. 
	Functions for shaping and manipulating arrays.
----------------------------------------------------------*/
//New implementation: An object that handles everything


//Malloc utility functions-----------------------
void *dan_malloc(int bytes, char const * name_var, char const* name_function){
	void *ret;
	if ((ret = malloc(bytes)) == NULL){
		//DAN_ERR_EXIT("Error allocating %s from function %s\n", name_var, name_function);
		throw std::bad_alloc();
	}
	return(ret);
}
void *dan_malloc(int bytes){
	void *ret;
	if ((ret = malloc(bytes)) == NULL){
		throw std::bad_alloc();
	}
	return(ret);
}

void dan_free(void* ptr){
	free(ptr);
}

//Creation of multidimensional arrays--------------

void *dan_flat2arrayND_ln(void* data, int size_elem, int dims, int *lengths){
	int size[20]; //max 20 dimensions (!). More would be too much, I think...
	char *indexes, *current, *next;
	int i, j, aux;
	
	if (dims == 1) return data;
	//compute the length of each level of indexes and its sum
	size[0] = lengths[0];
	aux = lengths[0];
	for (i = 1; i < dims -1; i ++){//is dims-1
		size[i] = size[i-1]*lengths[i];
		aux += size[i];
	}
	//allocate memory for the array of indexes. The last level (data) is already allocated!!
	indexes = (char*)dan_malloc(sizeof(void*)*aux, "indexes", "dan_flat2arrayND_ln");
	current = indexes; //current level to fill (leftmost index)
	
	//fill all levels except the one pointing directly to the data
	for (i = 0; i < dims -2; i++){
		next = current + size[i]*sizeof(void*); //next points to the next level of indexes
		//fill each element of the current level of indexes with the addresses
		for (j = 0; j < size[i]; j++){
			// of every lengths[i+1] elements of the next level.
			(*( (void**)current + j)) = next + lengths[i + 1] * j * sizeof(void*);
		}
		current = next;
	}
	//fill the last level (pointers to the real data!)
	for (j = 0; j < size[dims -2]; j++){
		*( (void**)current + j) = (char*)data + lengths[dims -1] * j * size_elem;
	}
	return (void*)indexes;
}

void *dan_flat2arrayND(void* data, int size_elem, int dims, ...){
	//read the variable input...
	int lengths[20];
	int i;

	std::va_list args;
	va_start(args,dims);
	for (i = 0; i < dims; i++){
		lengths[i] = va_arg(args, int);
	}
	va_end(args);
	return dan_flat2arrayND_ln(data, size_elem, dims, lengths);
}

void dan_free_arrayND_overhead(void* arreglo){
	free(arreglo); // just this!
}


//Versions for C++. These work with "new" instead of "malloc". Have to dealocate with "delete"
void *dan_flat2arrayND_cpp(void* data, int size_elem, const std::vector<int>& lengths){
	int dims = lengths.size();
	std::vector<int> size(dims);
	//void* indexes;
	unsigned char *current, *next;
	int i, j, aux;
	
	if (dims == 1) return data;
	//compute the length of each level of indexes and its sum
	size[0] = lengths[0];
	aux = lengths[0];
	for (i = 1; i < dims -1; i ++){//is dims-1
		size[i] = size[i-1]*lengths[i];
		aux += size[i];
	}
	//allocate memory for the array of indexes. The last level (data) is already allocated!!
	//indexes = (unsigned char*) ::operator new (sizeof(void*) * aux); //throw bad_alloc if failed.
	void* indexes =  ::operator new (sizeof(void*) * aux); //throw bad_alloc if failed.
	current = (unsigned char*) indexes; //current level to fill (leftmost index)
	
	//fill all levels except the one pointing directly to the data
	for (i = 0; i < dims - 2; i++){
		next = current + size[i] * sizeof(void*); //next points to the next level of indexes
		//fill each element of the current level of indexes with the addresses
		for (j = 0; j < size[i]; j++){
			// of every lengths[i+1] elements of the next level.
			(*( (void**)current + j)) = next + lengths[i + 1] * j * sizeof(void*);
		}
		current = next;
	}
	//fill the last level (pointers to the real data!)
	for (j = 0; j < size[dims -2]; j++){
		*( (void**)current + j) = (unsigned char*)data + lengths[dims -1]*j*size_elem;
	}
	//return (void*)indexes;
	return indexes;
}

void *dan_flat2arrayND_cpp(void* data, int size_elem, int dims, ...){
	//read the variable input...
	std::vector<int> lengths(dims);
	int i;
	std::va_list args;
	va_start(args,dims);
	for (i = 0; i < dims; i++){
		lengths[i] = va_arg(args, int);
	}
	va_end(args);
	return dan_flat2arrayND_cpp(data, size_elem, lengths);
}

/*Reimplementation of old functions*/
void* dan_flat2array2D(void* data, int size_elem, int dim_i, int dim_j){
	return dan_flat2arrayND(data, size_elem, 2, dim_i, dim_j);
}
void dan_free_array2D_overhead(void* myArray){
	dan_free_arrayND_overhead(myArray);
}

void* dan_flat2array3D(void* data, int size_elem, int dim_i, int dim_j, int dim_k){
	return dan_flat2arrayND(data, size_elem, 3, dim_i, dim_j, dim_k);
}
void dan_free_array3D_overhead(void* myArray,int dim_j){
	dan_free_arrayND_overhead(myArray);
}

//Contingency tables---------------

void dan_full_contingency_table(int **dataIJ, int &N, int &J, int* output, int *levelsJ, bool C_style, bool clear_output){
	//computes the contingency table
	//result in "output"
	//THIS IS A FULL EXPANDED CONTINGENCY TABLE. IT'S NOT A SPARSE REPRESENTATION. USE WITH CARE.
	if(clear_output) {
		//int len = dan_prod_vec(levelsJ, J);
		int len = 1;
		for (int j = 0; j < J; ++j) len *= *(levelsJ + J);
		memset(output, 0, sizeof(int)*len);
	}
	if(C_style){
		for (int i = 0; i < N; i++){
			output[ dan_index_2_linear_position_cstyle(dataIJ[i], J, levelsJ) ]++;
		}
	} else {
		for (int i = 0; i < N; i++){
			output[ dan_index_2_linear_position_fortranstyle(dataIJ[i], J, levelsJ) ]++;
		}
	}
}


