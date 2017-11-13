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

// Undef this when Cython no longer uses the deprecated API
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>     // This wants to be first on a Mac.

#include <everytrace.hpp>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/cython.hpp>
#include <exception>

// Cython exception handling
// http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#exceptions

namespace ibmisc {
namespace cython {

void init()
{
    everytrace_init();    // In case it wasn't initialized elsewhere; catch signals
    everytrace_exit = &everytrace_exit_exception;

    // Using Fortran dumps is OK
    // everytrace_dump_with_fortran = 0;
}


/** Makes sure a Numpy Array has the given type and dimension.  Raises
a Python exception if it does not.
@param type_num See http://docs.scipy.org/doc/numpy/reference/c-api.dtype.html eg: NPY_DOUBLE */
PyArrayObject *check_dimensions(
PyObject *ovec,
std::string const &vname,
int type_num,
int const *dims,
int ndim)
{
    // Check that it's not null
    if (!ovec) {
        // http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#exceptions
        (*ibmisc_error)(-1,
            "check_dimensions: Array object %s is null", vname.c_str());
    }


    // Check that it's type PyArrayObject
    if (!PyArray_Check(ovec)) {
// Alternately, we can throw our own exception.  This is about the same as using ibmisc_error.
//      throw std::invalid_argument(ibmisc::sprintf(
//          "check_dimensions: Object %s is not a Numpy array", vname.c_str()));
        (*ibmisc_error)(-1,
            "check_dimensions: Object %s is not a Numpy array", vname.c_str());
    }

    PyArrayObject *vec = (PyArrayObject *)ovec;

    // Check the data type and number of dimensions
    if (PyArray_DESCR(vec)->type_num != type_num || PyArray_NDIM(vec) != ndim)  {
        (*ibmisc_error)(-1,
            "check_dimensions: %s must be of type_num %d and %d dimensions (its is of type_num=%d and %d dimensions).", vname.c_str(), type_num, ndim, PyArray_DESCR(vec)->type_num, PyArray_NDIM(vec));
    }
//PyArrayObject *vec = (PyArrayObject *)ovec;

    // Check the dimensions themselves
    for (int i=0; i<ndim; ++i) {
        if (dims[i] < 0) continue;      // Don't check this dimension
        if (dims[i] != PyArray_DIM(vec, i)) {
            (*ibmisc_error)(-1,
                "%s: Array dimension #%d is %d, should be %d",
                vname.c_str(), i, PyArray_DIM(vec, i), dims[i]);
            // PyErr_SetString(PyExc_ValueError, buf);
        }
    }
    return vec;
}

}}  // namespace
