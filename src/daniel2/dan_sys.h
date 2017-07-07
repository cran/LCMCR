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

#ifndef _DAN_SYS_H
#define _DAN_SYS_H

#include <stdlib.h>
#include <stdio.h>
#include <cstdarg> //for ellipsis arguments
#ifdef _WIN32
#include <direct.h>
#define mkdir(x,y) _mkdir(x)
#else
#include <sys/stat.h>
#endif
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#endif

#ifndef __sun
#define restrict __restrict //Will work wigh gcc, g++, and VC++. Check this to make it compatible with more compilers
#else
#define restrict 
#endif



//Console output functions. With versions for R.
#ifdef USING_R
#define DAN_ERR_EXIT Rf_error
#define DAN_ERR_NOABORT Rf_warning
#define DAN_PRINTF Rprintf
#else
inline void DAN_ERR_EXIT(const char* f, ...) {
    std::va_list args;
    va_start(args,f);
    vfprintf(stderr, f,args);
	va_end(args);
	exit(1);
}
inline void DAN_ERR_NOABORT(const char* f, ...) {
    std::va_list args;
    va_start(args,f);
    //vfprintf(stderr, f, args);
	vprintf(f, args);
	va_end(args);
}
inline void DAN_PRINTF(const char* f, ...) {
    std::va_list args;
    va_start(args, f);
    vprintf(f, args);
	va_end(args);
}
#endif

#endif
