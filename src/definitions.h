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

#ifndef DEFINITIONS_H
#define	DEFINITIONS_H

/*Global definitions*/
#define GLOBAL_MAX_K 200 /*maximum K assumed for static allocations.*/
#define GLOBAL_MAX_EXP 70 /*maximum value for the argument of exp()*/

#ifdef _WIN32
#include <direct.h>
#define mkdir(x,y) _mkdir(x)
#else
#include <sys/stat.h>
#endif

#endif
