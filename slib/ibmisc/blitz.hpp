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

#ifndef IBMISC_BLITZ_HPP
#define IBMISC_BLITZ_HPP

#include <vector>
#include <array>
#include <blitz/array.h>
#include <blitz/tinyvec2.h>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/error.hpp>
#include <ibmisc/memory.hpp>

namespace ibmisc {

/** @defgroup blitz blitz.hpp
@brief Utilities for working with Blitz++ library.

@{
*/

template<class TypeT, int RANK>
bool is_row_major(blitz::Array<TypeT, RANK> &arr);

template<class TypeT, int RANK>
bool is_row_major(blitz::Array<TypeT, RANK> &arr)
{
    const blitz::TinyVector<int, RANK> &ordering(arr.ordering());

    for (int i=0; i<RANK; ++i) {
        if (ordering[RANK-i-1] != i) return false;
    }
    return true;
}

template<class TypeT, int RANK>
bool is_column_major(blitz::Array<TypeT, RANK> &arr);

template<class TypeT, int RANK>
bool is_column_major(blitz::Array<TypeT, RANK> &arr)
{
    const blitz::TinyVector<int, RANK> &ordering(arr.ordering());
    for (int i=0; i<RANK; ++i) {
        if (ordering[i] != i) return false;
    }
    return true;
}


/** Allocates an array in which all elements have the same value.
This is done by allocating a single value and setting all strides to 0.
By default, the array will have base of zero.  If a base of 1 is
desired, set storage = fortranArray. */
template<class TypeT, int RANK>
blitz::Array<TypeT, RANK> const_array(
    blitz::TinyVector<int,RANK> const &shape,
    TypeT const &val,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

template<class TypeT, int RANK>
blitz::Array<TypeT, RANK> const_array(
    blitz::TinyVector<int,RANK> const &shape,
    TypeT const &val,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    blitz::TinyVector<int,RANK> strides;
    for (int i=0;i<RANK;++i)
    {
        strides[i] = 0;
    }

    TypeT *data = new TypeT[1];
    *data = val;

    auto ret(blitz::Array<TypeT,RANK>(data,shape,strides,
        blitz::deleteDataWhenDone, storage));

    return ret;
}

/** Test if this array has all strides of 0.  This is useful because
Array.isStorageContiguous() returns false for such arrays. */
template<class TypeT, int RANK>
bool is_const_array(blitz::Array<TypeT, RANK> &arr)
{
    for (int i=0; i<RANK; ++i) if (arr.stride(i) != 0) return false;
    return true;
}
// ------------------------------------------------------------------
#define VECTOR_TO_BLITZ_BODY \
    blitz::TinyVector<int,1> shape(0); \
    blitz::TinyVector<int,1> strides(0); \
 \
    shape[0] = vec.size(); \
    strides[0] = 1;     /* Blitz++ strides in sizeof(T) units */ \
 \
    /* const_cast because Blitz++ can't construct a const Array */ \
    T *vecp = const_cast<T *>(&vec[0]); \
    return blitz::Array<T,1>(vecp, shape, strides, \
        blitz::neverDeleteData)


/** Converts a const std::vector to a const Blitz++ 1-D array that shares the same memory. */
template<class T>
blitz::Array<T,1> const to_blitz(std::vector<T> const &vec)
    { VECTOR_TO_BLITZ_BODY; }

template<class T>
blitz::Array<T,1> to_blitz(std::vector<T> &vec)
    { VECTOR_TO_BLITZ_BODY; }

/** Converts a const std::array to a const Blitz++ 1-D array that shares the same memory. */
#if 1    // C++14 only
template<class T, size_t LEN>
blitz::Array<T,1> const to_blitz(std::array<T, LEN> const &vec)
    { VECTOR_TO_BLITZ_BODY; }
#endif

template<class T, size_t LEN>
blitz::Array<T,1> to_blitz(std::array<T, LEN> &vec)
    { VECTOR_TO_BLITZ_BODY; }

#undef VECTOR_TO_BLITZ_BODY
// ------------------------------------------------------------------
template<class T>
std::vector<T> to_vector(blitz::Array<T,1> const &arr)
{
    std::vector<T> ret;
    for (size_t i=0; i < arr.shape()[0]; ++i) {
        ret.push_back(arr(i));
    }
    return ret;
}

template<class T, size_t len>
std::vector<T> to_vector(std::array<T,len> const &arr)
{
    std::vector<T> ret;
    for (size_t i=0; i < len; ++i) {
        ret.push_back(arr[i]);
    }
    return ret;
}

/** Cast data type (eg: float to double) while converting
from std::array to std::vector */
template<class SrcT, class DestT, size_t LEN>
std::vector<DestT> to_vector_cast(std::array<SrcT,LEN> const &arr)
{
    std::vector<DestT> ret;
    for (size_t i=0; i < LEN; ++i) {
        ret.push_back(arr[i]);
    }
    return ret;
}


// ------------------------------------------------------------------
template<class TinyT, class ArrayT, int RANK>
blitz::TinyVector<TinyT, RANK> to_tiny(std::array<ArrayT, RANK> const &arr)
{
    blitz::TinyVector<TinyT, RANK> ret;
    for (int k=0; k<RANK; ++k) ret[k] = arr[k];
    return ret;
}

template<class TinyT, class ArrayT, int RANK>
void to_tiny(blitz::TinyVector<TinyT, RANK> &ret, std::array<ArrayT, RANK> const &arr)
{
    for (int k=0; k<RANK; ++k) ret[k] = arr[k];
}
// ------------------------------------------------------------------


template<class ArrayT, class TinyT, int RANK>
std::array<ArrayT, RANK> to_array(blitz::TinyVector<TinyT, RANK> const &tiny)
{
    std::array<ArrayT, RANK> ret;
    for (int k=0; k<RANK; ++k) ret[k] = tiny[k];
    return ret;
}

/** Convert std::vector to std::array */
template<class ArrayT, class VecT, int RANK>
std::array<ArrayT, RANK> to_array(std::vector<VecT> const &vec)
{
    std::array<ArrayT, RANK> ret;
    for (int k=0; k<RANK; ++k) ret[k] = vec[k];
    return ret;
}

template<class ArrayT, class TinyT, int RANK>
void to_array(std::array<ArrayT, RANK> &ret, blitz::TinyVector<TinyT, RANK> const &tiny)
{
    for (int k=0; k<RANK; ++k) ret[k] = tiny[k];
}


// ------------------------------------------------------------------

// =============================================================
/** Frees memory associated with a blitz::Array */
template<class T, int len>
void free_array(blitz::Array<T, len> &array)
{
    array.reference(blitz::Array<T,len>(0, blitz::shape(0), blitz::neverDeleteData));
}


/** Makes sure a blitz::Array dimension.
Raises a Python exception if it does not. */

/** Checks that the dimensions of an array are what we think they
should be.  This is used in Python interface code to add sophisticated
type checking to functions.

@param vname Name of this variables (used in error comments)
@param arr The array to check dimensions
@param dims The expected dimensions.  If a dims[i] < 0, then dimension i is not checked. */
template<class T, int rank>
void check_dimensions(
std::string const &vname,
blitz::Array<T, rank> const &arr,
std::vector<int> const &dims)
{
    for (int i=0; i<rank; ++i) {
        if (dims[i] >= 0 && arr.extent(i) != dims[i]) {
            (*ibmisc::ibmisc_error)(-1,
                "Error in %s: expected dimension #%d = %d (is %d instead)\n",
                vname.c_str(), i, dims[i], arr.extent(i));
        }
    }
}
// ------------------------------------------------------------

/** Changes a C-style Blitz++ array (biggest stride in first dimension
and zero-based indexing) to a Fortran-style Blitz++ array (biggest
stride in last dimension and one-based indexing) that shares the same memory. */
template<class T, int rank>
blitz::Array<T, rank> c_to_f(blitz::Array<T, rank> &arr);

template<class T, int rank>
blitz::Array<T, rank> c_to_f(blitz::Array<T, rank> &arr)
{
    // Initialize an 11-dim vector of transpositions
    // (because transpose() doesn't take a TinyVector)
    int const max_dims = 11;
    int rev[max_dims];
    for (int i=rank; i<max_dims; ++i) rev[i] = 0;

    // Reverse dimensions
    for (int i=0; i<rank; ++i) rev[i] = rank-i-1;
    auto ret(arr.transpose(rev[0], rev[1], rev[2], rev[3], rev[4], rev[5], rev[6], rev[7], rev[8], rev[9], rev[10]));

    // Re-base to 1
    blitz::TinyVector<int, rank> base(1);
    ret.reindexSelf(base);

    return ret;
}
// ------------------------------------------------------------
/** Changes a C-style Blitz++ array (biggest stride in first dimension
and zero-based indexing) to a Fortran-style Blitz++ array (biggest
stride in last dimension and one-based indexing) that shares the same memory. */
template<class T, int rank>
blitz::Array<T, rank> f_to_c(blitz::Array<T, rank> &arr);

template<class TypeT, int RANK>
inline blitz::Array<TypeT, RANK> f_to_c(blitz::Array<TypeT, RANK> &arr)
{
    // Set up shape and strides
    blitz::TinyVector<int,RANK> shape;
    blitz::TinyVector<int,RANK> strides;
    blitz::GeneralArrayStorage<RANK> storage;
    for (int i=0;i<RANK;++i)
    {
        int j = RANK-1 - i;
        shape[i] = arr.extent(j);
        // Python/Numpy strides are in bytes, Blitz++ in sizeof(T) units.
        strides[i] = arr.stride(j);
        storage.base()[i] = arr.lbound(j)-1;
        // Ordering is not needed because we're using stride
        // storage.ordering()[i] = i;      // Fortran ordering, blitz++ manual p31
    }

    auto ret(blitz::Array<TypeT,RANK>(arr.data(),shape,strides,
        blitz::neverDeleteData, storage));

    return ret;

}

template<class TypeT, int RANK>
blitz::Array<TypeT, RANK> const f_to_c(blitz::Array<TypeT, RANK> const &arr)
{
    return f_to_c(
        *const_cast<blitz::Array<TypeT, RANK> *>(&arr));
}
// ---------------------------------------------------------
/** Allows negative dimensions in dest_shape, as with Python's numpy.reshape().
See https://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.reshape.html#numpy.reshape */

/** Helper function to determine array shapes for reshape() below. */
template<int src_ndim, int dest_ndim>
blitz::TinyVector<int,dest_ndim> reshape_get_dest_shape(
    blitz::TinyVector<int,src_ndim> const &src_shape,
    blitz::TinyVector<int,dest_ndim> const &dest_shape);

template<int src_ndim, int dest_ndim>
blitz::TinyVector<int,dest_ndim> reshape_get_dest_shape(
    blitz::TinyVector<int,src_ndim> const &src_shape,
    blitz::TinyVector<int,dest_ndim> const &dest_shape)
{
    blitz::TinyVector<int,dest_ndim> ret;

    /* Check dimensions */
    long src_n = 1;
    for (int i=0; i<src_ndim; ++i) {
        src_n *= src_shape[i];
    }

    int negi = -1;
    long dest_n = 1;
    for (int i=0; i<dest_ndim; ++i) {
        if (dest_shape[i] < 0) {
            if (negi >= 0) (*ibmisc::ibmisc_error)(-1,
                "reshape() destination cannot have more than one dimension <0 (but it has at least dimensions #%d and #%d)", negi, i);
            negi = i;
        } else {
            dest_n *= dest_shape[i];
            ret[i] = dest_shape[i];
        }
    }

    // Fill in missing dimension
    if (negi >= 0) {
        ret[negi] = src_n / dest_n;    // integer divide
        if (ret[negi] * dest_n != src_n) (*ibmisc::ibmisc_error)(-1,
            "Non-integral missing dimension #i: %d / %d", negi, src_n, dest_n);
    } else if (src_n != dest_n) {
        (*ibmisc::ibmisc_error)(-1,
            "blitz.hpp, ibmisc::reshape(): Total dimension mismatch, src=%ld, dest=%ld\n", src_n, dest_n);
    }

    return ret;
}


/** Reshape an array, to another one with different dimensions but the
    same number of elements.  Returns an array aliasing the same
    internal memory.

    To reshape `blitz::Array<double,2> A` to rank 1:
       reshape<double,2,1>(A, {-1});

@param T Array type
@param src_ndim Rank of the array to be reshaped
@param dest_ndim Rank of the destination array to return
@param src The array to be reshaped
@param dest_shape Shape of the array to return.
    The total size implied by this shape (all elements multipled
    together) must equal the size of src.  Up to one dimension may be
    -1: as with Python's numpy.reshape(), that dimension will be set
    to "all remaining values."
@return The reshaped array. */
template<class T, int src_ndim, int dest_ndim>
extern blitz::Array<T, dest_ndim> reshape(
    blitz::Array<T, src_ndim> &src,
    blitz::TinyVector<int,dest_ndim> const &dest_shape)
{
    /* Do the reshaping */
    return blitz::Array<T,dest_ndim>(src.data(),
        reshape_get_dest_shape(src.shape(), dest_shape),
        blitz::neverDeleteData);
}

template<class T, int src_ndim, int dest_ndim>
extern blitz::Array<T, dest_ndim> const reshape(
    blitz::Array<T, src_ndim> const &src,
    blitz::TinyVector<int,dest_ndim> const &dest_shape)
{
    /* Do the reshaping */
    return blitz::Array<T,dest_ndim>(const_cast<T *>(src.data()),
        reshape_get_dest_shape(src.shape(), dest_shape),
        blitz::neverDeleteData);
}

/** Simple method to reshape to a 1-D array.
@param storage Set to blitz::fortranArray to get 1-based indexing.
@param arr The array to reshape
@param lbound Lower bound to use on reshaped array: 0 (C) or 1 (Fortran) */
template<class TypeT, int RANK>
blitz::Array<TypeT,1> reshape1(
    blitz::Array<TypeT, RANK> &arr,
    int lbound = 0);

template<class TypeT, int RANK>
blitz::Array<TypeT,1> reshape1(blitz::Array<TypeT, RANK> &arr,
    int lbound)
{
    if (!(is_const_array(arr) || arr.isStorageContiguous())) (*ibmisc_error)(-1,
        "Array must be contiguous for reshape1().");

    // Get the stride we will use
    blitz::diffType min_stride = arr.stride(0);
    for (int i=1; i<RANK; ++i) min_stride = std::min(min_stride, arr.stride(i));

    blitz::TinyVector<int,1> shape, stride;
    blitz::GeneralArrayStorage<1> stor;
    shape[0] = arr.size();
    stride[0] = min_stride;
    stor.base()[0] = lbound;
    // stor Ordering is not needed because we're using stride
    return blitz::Array<TypeT,1>(arr.data(), shape, stride,
        blitz::neverDeleteData, stor);
}

/** const version of reshape1() */
template<class TypeT, int RANK>
blitz::Array<TypeT,1> const reshape1(blitz::Array<TypeT, RANK> const &arr,
    int lbound = 0)
{
    return reshape1(
        *const_cast<blitz::Array<TypeT, RANK> *>(&arr),
        lbound);
}

// These templates SHOULD work.  But they haven't been tested or used,
// so they're commented out for now.
// ------------------------------------------------
#if 0
template<class T, int len>
blitz::TinyVector<T, len> vector_to_tiny(std::vector<T> const &vec)
{
    if (vec.size() != len) {
        (*ibmisc::ibmisc_error)(-1,
            "vector_to_tiny(): vector length %ld does not match declared length %d\n", vec.size(), len);
    }

    blitz::TinyVector<T, len> ret;
    for (int i=0; i < len; ++i) {
        ret[i] = vec[i];
    }
    return ret;
}
#endif
// ------------------------------------------------
template<class T, int len>
std::vector<T> tiny_to_vector(blitz::TinyVector<T, len> const &tiny)
{
    std::vector<T> ret;
    ret.reserve(len);
    for (int i=0; i<len; ++i) ret.push_back(tiny[i]);
    return ret;
}
// ----------------------------------------------------------------
#if 0
/** Casts from one Blitz array type to another */
template<class SrcT, class DestT, int RANK>
blitz::Array<DestT, RANK> blitz_cast(blitz::Array<SrcT, RANK> const &src)
{
    if (!src.isStorageContiguous) (*ibmisc_error)(-1,
        "Storage must be contiguous");

    blitz::Array<DestT, RANK> dest(src.shape());

    auto srci(src.begin());
    auto desti(dest.begin());
    for (; srci != src.end(); ++srci, ++desti) {
        *desti = (DestT)*srci;
    }
    return dest;
}

template<class SrcT, int RANK>
blitz::Array<double, RANK> to_double(blitz::Array<SrcT, RANK> const &src)
    { return blitz_cast<SrcT, double, RANK>(src); }

#endif
// ----------------------------------------------------------------

/** @} */
}   // NAMESPACE

#endif
