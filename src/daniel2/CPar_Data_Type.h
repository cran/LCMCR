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

#ifndef _CPar_Data_Type_H_
#define _CPar_Data_Type_H_

#include<string>

//A class for descrribing data types. 
//It helps to keep track of data characteristics when performing untyped operations.

class CPar_Data_Type {
public:
	typedef enum {T_INT, T_DOUBLE, T_BYTE} data_type_t;
	CPar_Data_Type(data_type_t type) : type(type){
		switch (type){
		case T_INT:
			name = "int";
			n_bytes = sizeof(int);
			break;
		case T_DOUBLE:
			name = "double";
			n_bytes = sizeof(double);
			break;
		case T_BYTE:
			name = "byte";
			n_bytes = sizeof(unsigned char);
			break;
		}
	}
	//use default copy constructor and "=" operator.
	int get_bytes(){ return n_bytes; }
	std::string& get_type_name(){ return name; }
	data_type_t get_data_type() { return type; }
private:
	int n_bytes;
	std::string name;
	data_type_t type;
};

#endif
