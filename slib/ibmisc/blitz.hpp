#ifndef SPSPARSE_BLITZ_HPP
#define SPSPARSE_BLITZ_HPP

#include <vector>
#include <blitz/array.h>

namespace ibmisc {

/** @defgroup blitz blitz.hpp
@brief Utilities for working with Blitz++ library.

@{
*/

// ------------------------------------------------------------------
#define VECTOR_TO_BLITZ_BODY \
    blitz::TinyVector<int,1> shape(0); \
    blitz::TinyVector<int,1> strides(0); \
 \
	shape[0] = vec.size(); \
	strides[0] = 1;		/* Blitz++ strides in sizeof(T) units */ \
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

/** @} */
}	// NAMESPACE

#endif
