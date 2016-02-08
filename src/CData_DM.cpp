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

/*
	CData_DM class implementation
	2010-10-20 - multilevel extensions
	2011-06-02 - Contingency table generation
	2013-10-16 - Redesign. Using CParams_generic container and new CData base class.
*/

#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include "CData_DM.h"

#include "daniel2/dan_array_utils.h"
#include "daniel2/dan_math.h"
//using namespace std;







