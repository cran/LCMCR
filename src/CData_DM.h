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

#if !defined(_CDATA_DM_H)
#define _CDATA_DM_H

#include "daniel2/CData.h"
#include "daniel2/CPar_Data_Type.h"
//#include <iostream>
class CData_DM : public CData{
public:
	CData_DM(){
		Declare();
		ncells = 0;
	}
	void Set_Manually(int *x_flat, int J, int n, int *levels){
		_Load_Variable("x", x_flat, n, J);
		_Load_Variable("levelsJ", levels, J);
	}
	//data storage. These should be private, but this is going to be more efficient.
	int ** x; //responses vector
	int *levelsJ; // (levelsJ[j] = #levels of variable j)
	int L; //max #of levels (for dimentions)
	int J;
	int n; //sample size
	int **cells; //Occupied cells in the contingency table (Initialized in TabulateContingency())
	int *freqM; //Frequencies of counts in contingency tabl (Initialized in TabulateContingency()).
	int ncells; //number of occupied cells in the contingency table (Initialized in TabulateContingency())
protected:
	void Close_Loading(){
		CData::Close_Loading();
		_set_internal_pointer("x", &x);
		_set_internal_pointer("levelsJ", &levelsJ);
		//compute statistics
		J = data_container["x"].get_dim_lengths()[1];
		n = data_container["x"].get_dim_lengths()[0];
		L = *std::max_element(levelsJ, levelsJ + J);
		ncells = 0;//THe contingency table is only computed if requested explicitly using TabulateContingency()
		TabulateContingency();
	}
private:
	void Declare(){
		_Declare_Variable("x", CPar_Data_Type::T_INT, 2, false);
		_Declare_Variable("levelsJ", CPar_Data_Type::T_INT, 1, false);
	}
public:
	void TabulateContingency(){
		if (!is_loaded()){
			throw std::runtime_error("Cannot compute contingency table. Data not read.");
		}
		typedef std::map<std::vector<int>, int> tabla;
		tabla s;
		//tabulate the data (cells in tabla.first, counts in tabla.second)
		for (int i = 0; i < n; i++){// int* it = x[0]; it != x[0] + n * J; it += J){
			++s[std::vector<int>(x[i], x[i] + J)];
		}
		ncells = s.size();
		//allocate space
		cells = (int**)_Declare_and_Allocate_derived("cells", CPar_Data_Type::T_INT, 2, ncells, J);
		freqM = (int*)_Declare_and_Allocate_derived("freqM", CPar_Data_Type::T_INT, 1, ncells);
		//and copy the result.
		int m = 0;
		for (tabla::iterator it = s.begin(); it != s.end(); ++it, ++m){
			std::copy(&(it->first[0]), &(it->first[0]) + J, cells[m]);
			freqM[m] = it->second;
		}
	}
	bool Contingency_Tabulated(){ return ncells > 0; }
};

#endif  //_CDATA_DM_H
