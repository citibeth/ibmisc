/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IBMISC_ERROR_HPP
#define IBMISC_ERROR_HPP

#include <everytrace.hpp>

/** @defgroup ibmisc ibmisc.hpp
@brief Basic stuff common to all ibmisc */
namespace ibmisc {

typedef everytrace::Exception Exception;

/** Use Everytrace error handler & back-end by default; user or other
    library can change if needed. */
extern everytrace_error_ptr ibmisc_error;

}   // namespace
/** @} */

#endif // Guard
