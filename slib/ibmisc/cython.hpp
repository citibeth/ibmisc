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
// ======================================================

// Primitive Python types
class pytype_long {};

// ======================================================

/** General template to box a primitive value into a Python object */
template<class PyTypeT, class TypeT>
inline PyObject *py_build(TypeT const &ival)
{
    PyErr_SetString(PyExc_ValueError, "py_build(): Unknown type combination");
    ibmisc::cython::exit();
}
template<> inline PyObject *py_build<pytype_long, int>(int const &ival)
    { return PyLong_FromLong(ival); }
template<> inline PyObject *py_build<pytype_long, long>(long const &ival)
    { return PyLong_FromLong(ival); }


// ======================================================

/** Convenient template to build a Python tuple from a std::array
http://effbot.org/pyfaq/how-do-i-use-py-buildvalue-to-create-a-tuple-of-arbitrary-length.htm */
template<class PyTypeT, class TypeT, int RANK>
PyObject *py_build_tuple(std::array<TypeT,RANK> val);

template<class PyTypeT, class TypeT, int RANK>
PyObject *py_build_tuple(std::array<TypeT,RANK> val)
{
    PyObject *result = PyTuple_New(RANK);
    if (!result) return NULL;

    for (int i=0; i<RANK; ++i) {
        PyObject *value = py_build<PyTypeT,TypeT>(val[i]);
        if (!value) {
            Py_DECREF(result);
            return NULL;
        }
        PyTuple_SetItem(result, i, value);
    }
    return result;
}
// ======================================================

// --------------------------------------------------------------------
// Convert C++ template types to Numpy type_nums

template<class T>
inline int np_type_num()
{
    PyErr_SetString(PyExc_ValueError, "np_type_num(): Unknown type_num");
    ibmisc::cython::exit();
}

template<> inline int np_type_num<char>()
    { return NPY_INT8; }

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
@param dims Dimensions are checked to be this size.  Use -1 to accept whatever
            dimension the Numpy array comes with.
*/
template<class T, int N>
blitz::Array<T,N> np_to_blitz(
PyObject *ovec,
std::string const &vname,
std::array<int,N> dims,
blitz::GeneralArrayStorage<N> const &storage = blitz::GeneralArrayStorage<N>());

template<class T, int N>
blitz::Array<T,N> np_to_blitz(
PyObject *ovec,
std::string const &vname,
std::array<int,N> dims,
blitz::GeneralArrayStorage<N> const &storage)
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
        blitz::neverDeleteData, storage);
}


/** Allocates a Numpy array of a given type and size
For memory allocation issues, see:
http://grokbase.com/t/python/python-list/005ttee301/decrefing-and-pyarray-functions-in-c-extensions

    Every PyArrayObject that you successfully receive from a PyArray_XXXXXX
    call has an increased reference count (you own a reference). Unless you
    are returning the object to the user (using return PyArray_Return(obj) or
    PyArray_BuildValue with the "N" construct) you must DECREF the object
    before leaving the subroutine or you will create a memory leak..
*/
template<class T, int N>
PyObject *new_pyarray(std::array<int,N> dims)
{
    return PyArray_FromDims(N, &dims[0], np_type_num<T>());
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
PyObject *copy_blitz_to_np(blitz::Array<T,N> const &array)
{
    std::array<int,N> dims;
    for (int i=0; i<N; ++i) dims[i] = array.extent(i);
    PyObject *array_py = PyArray_FromDims(N, &dims[0], np_type_num<T>());

    // Wrap as blitz for convenience
    auto array_bl(np_to_blitz<T,N>(array_py, "array_py", dims));

    // Copy from the original Blitz array to the Blitzified Numpy array
    array_bl = array;
    return array_py;
}



template<class ArrayT>
PyObject *spsparse_to_tuple(ArrayT const &A);

template<class ArrayT>
PyObject *spsparse_to_tuple(ArrayT const &A)
{
    const int rank = ArrayT::rank;
    typedef typename ArrayT::index_type index_type;
    typedef typename ArrayT::val_type val_type;

    // Allocate Python arrays and Blitzify
    std::array<PyObject *,rank> indices_np;
    std::array<blitz::Array<index_type,1>, ArrayT::rank> indices_bl;
    std::array<int,1> dims {(int)A.size()};
    for (int k=0; k<ArrayT::rank; ++k) {
        indices_np[k] = PyArray_FromDims(
            1, &dims[0], np_type_num<index_type>());
        indices_bl[k].reference(np_to_blitz<index_type,1>(indices_np[k], "index", dims));
    }
    PyObject *values_np = PyArray_FromDims(
        1, &dims[0], np_type_num<typename ArrayT::val_type>());
    auto values_bl(np_to_blitz<val_type,1>(values_np, "values", dims));

    // Copy into the (Numpy) arrays we just allocated
    int i=0;
    for (auto ii=A.begin(); ii != A.end(); ++ii,++i) {
        for (int k=0; k<ArrayT::rank; ++k)
            indices_bl[k](i) = ii->index(k);
        values_bl(i) = ii->value();
    }

    // Bind it up in a tuple
    PyObject *indices_t = PyTuple_New(ArrayT::rank);
    PyObject *shape_t = PyTuple_New(ArrayT::rank);
    for (int k=0; k<ArrayT::rank; ++k) {
        // blitz::Array<typename ArrayT::val_type, ArrayT::rank> idx(indices(k));
        PyTuple_SetItem(indices_t, k, indices_np[k]);
        PyTuple_SetItem(shape_t, k, Py_BuildValue("i", A.shape()[k]));
    }


    PyObject *data_t = Py_BuildValue("OO", values_np, indices_t);

    // For scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))
    PyObject *ret = Py_BuildValue("OO", data_t, shape_t);

    return ret;
}


}}  // namespace
