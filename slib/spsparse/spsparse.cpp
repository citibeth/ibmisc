/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <functional>
#include <cstdlib>
#include <spsparse/spsparse.hpp>
#include <exception>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif


namespace spsparse {

void default_error(int retcode, const char *format, ...)
{
    va_list arglist;

    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);
    fprintf(stderr, "\n");

#ifdef USE_EVERYTRACE
    everytrace_exit(-1);
#endif
    throw spsparse::Exception();
//  exit(-1);
}

error_ptr spsparse_error = &default_error;

const std::array<int,2> ROW_MAJOR = {0,1};
const std::array<int,2> COL_MAJOR = {1,0};

}   // Namespace
