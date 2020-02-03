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

#ifndef _CDATA_DM_STRAT_H
#define _CDATA_DM_STRAT_H

#include"CData_DM.h"

class CData_DM_Strat : public CData_DM{
public:
	int C; //Number of strata
	int *nC; //Counts of observed ind by stratum
	int** cells_strat;
	int* freqM_strat;
	int ncells_strat;
	int** cells_mis;
	int* freqM_mis;
	int ncells_mis;

	CData_DM_Strat(){
		//Initialize
		ncells_strat = 0;
		ncells_mis = 0;
	}
protected:
	void Close_Loading(){
		CData_DM::Close_Loading();

		//Calculate Statistics.
		Process_Tables();
	}
private:
	void Declare(){
		//Additional declaration
	}
public:
	void Process_Tables(){
		typedef std::map<std::vector<int>, int> tabla;
		if (!is_loaded()){
			throw std::runtime_error("Cannot process data. Data not yet loaded.");
		}
		this->C = this->levelsJ[0] - 1; //Number of strata

		//tabulate the data (cells in tabla.first, counts in tabla.second)
		tabla s_strat, s_mis;
		std::vector<int> v(J);
		for (int i = 0; i < n; i++){
			if(x[i][0] == C){
				//If missing stratification
				std::copy(x[i], x[i] + J, &(v[0]));
				v[0] = -1; //Make sure to have an invalid value for missing data.
				++s_mis[ std::vector<int>(v.begin(), v.end()) ];
			} else {
				//If regular data
				++s_strat[ std::vector<int>(x[i], x[i] + J) ];
			}
		}
//		for (tabla::iterator it = s_mis.begin(); it != s_mis.end(); it++){
//			for (int i = 0; i < levelsJ[0]; i++){
//				std::copy(it->first.begin(), it->first.end(), v.begin());
//				v[0] = i;
//			}
//		}

		//allocate space
		ncells_strat = s_strat.size();
		ncells_mis = s_mis.size();
		nC = (int *)_Declare_and_Allocate_derived("nC", CPar_Data_Type::T_INT, 1, C);
		cells_strat = (int**)_Declare_and_Allocate_derived("cells_strat", CPar_Data_Type::T_INT, 2, ncells, J);
		freqM_strat = (int*)_Declare_and_Allocate_derived("freqM_strat", CPar_Data_Type::T_INT, 1, ncells);
		cells_mis = (int**)_Declare_and_Allocate_derived("cells_mis", CPar_Data_Type::T_INT, 2, ncells_mis, J);
		freqM_mis = (int*)_Declare_and_Allocate_derived("freqM_mis", CPar_Data_Type::T_INT, 1, ncells_mis);

		//and copy the results.
		int m = 0; 
		(*this)["nC"].fill(0);
		for (tabla::iterator it = s_strat.begin(); it != s_strat.end(); ++it, ++m){
			std::copy(&(it->first[0]), &(it->first[0]) + J, cells_strat[m]);
			freqM_strat[m] = it->second;
			int c = it->first[0];
			nC[c] += it->second;

		}
		m = 0;
		for (tabla::iterator it = s_mis.begin(); it != s_mis.end(); ++it, ++m){
			std::copy(&(it->first[0]), &(it->first[0]) + J, cells_mis[m]);
			freqM_mis[m] = it->second;
		}
		//Maybe index repeated elements. 
	}

	bool data_processed(){ return ncells > 0; }
};

#endif  //_CDATA_DM_H
