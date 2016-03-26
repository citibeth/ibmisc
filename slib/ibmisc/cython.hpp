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

#pragma once

#include <Python.h>
#include <numpy/arrayobject.h>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/ibmisc.hpp>


#include <blitz/array.h>

/** Help for Cython bindings to things in this library. */

namespace ibmisc {
namespace cython {

/** Call when initializing the Cython extension. */
void init();

extern PyArrayObject *check_dimensions(
PyObject *ovec,
std::string const &vname,
int type_num,
int const *dims,
int ndim);

// -------------------------------------------------

static inline ibmisc::NcIO *new_ncio(std::string fname, std::string sfileMode)
{
    netCDF::NcFile::FileMode fileMode;
    if (sfileMode == "read") fileMode = netCDF::NcFile::FileMode::read;
    else if (sfileMode == "write") fileMode = netCDF::NcFile::FileMode::write;
    else if (sfileMode == "replace") fileMode = netCDF::NcFile::FileMode::replace;
    else if (sfileMode == "newFile") fileMode = netCDF::NcFile::FileMode::newFile;
    else {
        (*ibmisc::ibmisc_error)(-1,
            "Bad file mode for netCDF::NcFile::FileMode: %s", sfileMode.c_str());
    }


    return new ibmisc::NcIO(fname, fileMode);
}

inline void exit()
{
    // Not sure what we should do here...
    // We want to exit the Cython code and return an exception to Python
}

// --------------------------------------------------------------------
// Convert template types to Numpy type_nums

template<class T>
inline int np_type_num()
{
    PyErr_SetString(PyExc_ValueError, "np_type_num(): Unknown type_num");
    ibmisc::cython::exit();
}

template<> inline int np_type_num<short>()
    { return NPY_SHORT; }

template<> inline int np_type_num<int>()
    { return NPY_INT; }

template<> inline int np_type_num<long>()
    { return NPY_LONG; }

template<> inline int np_type_num<float>()
    { return NPY_FLOAT; }

template<> inline int np_type_num<double>()
    { return NPY_DOUBLE; }

template<> inline int np_type_num<long double>()
    { return NPY_LONGDOUBLE; }

template<> inline int np_type_num<long long>()
    { return NPY_LONGLONG; }

/** TODO: Add the following basic Numpy types

#define NPY_SIZEOF_COMPLEX_FLOAT 8
#define NPY_SIZEOF_DOUBLE 8
#define NPY_SIZEOF_COMPLEX_DOUBLE 16
#define NPY_SIZEOF_LONGDOUBLE 16
#define NPY_SIZEOF_COMPLEX_LONGDOUBLE 32
#define NPY_SIZEOF_PY_INTPTR_T 8
#define NPY_SIZEOF_PY_LONG_LONG 8
#define NPY_SIZEOF_LONGLONG 8
*/
// --------------------------------------------------------------------

/** Convert a Numpy array to a blitz one, using the original's data (no copy).
N is the rank of the array.
@see: http://mail.scipy.org/pipermail/numpy-discussion/2004-March/002892.html

Use this template to convert from PyObject to a blitz::Array of a specific
type and dimensions.

Example:
    PyObject *Val_py;
    auto Val(np_to_blitz<double,3>(Val_py, "Val", {3, -1, 4}));

@param type_num Must match T
*/
template<class T, int N>
blitz::Array<T,N> np_to_blitz(
PyObject *ovec,
std::string const &vname,
std::array<int,N> dims);

template<class T, int N>
blitz::Array<T,N> np_to_blitz(
PyObject *ovec,
std::string const &vname,
std::array<int,N> dims)
{
    // Check data types and cast
    PyArrayObject *vec = check_dimensions(ovec, vname, np_type_num<T>(), &dims[0], N);

    int const T_size = sizeof(T);
    assert(T_size == PyArray_ITEMSIZE(vec));

    blitz::TinyVector<int,N> shape(0);
    blitz::TinyVector<int,N> strides(0);

    // Set up shape and strides
    for (int i=0;i<N;++i) {
        shape[i]   = PyArray_DIM(vec, i);
        // Python/Numpy strides are in bytes, Blitz++ in sizeof(T) units.
        strides[i] = PyArray_STRIDE(vec, i) / T_size;
    }

    return blitz::Array<T,N>((T*) PyArray_DATA(vec),shape,strides,
        blitz::neverDeleteData);
}





#if 0
template<class T, int N>
PyObject *blitz_to_np(
blitz::Array<T,N> const &arr,
std::string const &vname,
std::array<int,N> const &dims)
{
    /* See:
http://gael-varoquaux.info/programming/cython-example-of-exposing-c-computed-arrays-in-python-without-data-copies.html
http://docs.scipy.org/doc/numpy-1.10.0/reference/c-api.array.html
    */
}
#endif

template<class T, int N>
PyObject *blitz_to_np(blitz::Array<T,N> const &array)
{
    std::array<int,N> dims;
    for (int i=0; i<N; ++i) dims[i] = array.extent(i);
    PyObject *array_py = PyArray_FromDims(N, &dims[0], np_type_num<T>());

    // Wrap as blitz for convenience
    auto array_bl(np_to_blitz<T,N>(array_py, "array_py", dims));

    // Copy from the arraytor to the Blitz array
    array_bl = array;
    return array_py;
}



template<class ArrayT>
PyObject *spsparse_to_tuple(ArrayT A);

template<class ArrayT>
PyObject *spsparse_to_tuple(ArrayT A)
{
    PyObject *indices_t = PyTuple_New(ArrayT::rank);
    PyObject *shape_t = PyTuple_New(ArrayT::rank);
    for (int k=0; k<ArrayT::rank; ++k) {
        // blitz::Array<typename ArrayT::val_type, ArrayT::rank> idx(indices(k));
        PyTuple_SetItem(indices_t, k, blitz_to_np<typename ArrayT::index_type,1>(A.indices(k)));
        PyTuple_SetItem(shape_t, k, Py_BuildValue("i", A.shape[k]));
    }

    PyObject *data_t = Py_BuildValue("OO",
        blitz_to_np<typename ArrayT::val_type, 1>(A.vals()),
        indices_t);

    // For scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))
    PyObject *ret = Py_BuildValue("OO", data_t, shape_t);

    return ret;
}


}}  // namespace
