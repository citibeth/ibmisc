#ifndef SPSPARSE_ARRAY_HPP
#define SPSPARSE_ARRAY_HPP

#include <array>
#include <cstdio>
#include <functional>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <ibmisc/iter.hpp>
#include <ibmisc/blitz.hpp>

#include <spsparse/spsparse.hpp>
#include <spsparse/algorithm.hpp>
#include <spsparse/xiter.hpp>
#include <spsparse/accum.hpp>

namespace spsparse {

/** @defgroup array array.hpp
@brief Basic sparse arrays

@{
*/


// -----------------------------------------------------
/** @brief Select out just one dimension of the index on iteration.

Wraps VectorCooArray::iterator, changing operator*() to produce produces the
index in just one dimension.

Code Example
@code
VectorCooMatrix<int, double> const A;
typedef DimIndexIter<decltype(A)::const_iterator> DIType;
for (DIType ii(1, A.begin()); ii != DIType(1, A.end()); ++ii)
	printf("Element with column %d and value %d\n", *ii, ii.val());
@endcode

@see spsparse::VectorCooMatrix::dim_iter(), spsparse::VectorCooMatrix::dim_begin(), spsparse::VectorCooMatrix::dim_end()
*/
template<class ValueT, class ValT, class IterT>
class DimIndexIter : public ibmisc::forward_iterator<ValueT, DimIndexIter<ValueT, ValT, IterT>>
{
public:
	IterT wrapped;
	const int dim;

	DimIndexIter(int _dim, IterT const &&ii) : wrapped(ii), dim(_dim) {}

	ValueT operator*()
		{ return wrapped.index(dim); }
	ValT &val()
		{ return wrapped.val(); }
	ValT const val() const
		{ return wrapped.val(); }

	DimIndexIter &operator++()
		{ ++wrapped; return *this; }
	bool operator==(const DimIndexIter& rhs) const
		{return wrapped == rhs.wrapped;}
};
// -----------------------------------------------------
template<class IndicesT, class IterIndexT, int RANK, class IterValT, class CollectionT>
class CooIterator
{
protected:
	CollectionT * const parent;
	int i;
public:
	static const int rank = RANK;
	typedef IndicesT indices_type;
	typedef indices_type value_type;	// Standard STL: Type we get upon operator*()
	typedef IterIndexT index_type;	// Our convention
	typedef IterValT val_type;	// Extension: Type we get upon val()

	CooIterator(CollectionT *p, int _i) : parent(p), i(_i) {}


	indices_type operator[](int n)
		{ return parent->index(i+n); }
	indices_type index()
		{ return parent->index(i); }
	indices_type operator*()
		{ return this->operator[](0); }

	CooIterator &operator+=(int n)
		{i += n; return *this; }
	CooIterator& operator--()
		{ return this->operator+=(-1); }
	CooIterator& operator++()
		{ return this->operator+=(1); }
	CooIterator &operator-=(int n)
		{ return this->operator+=(-n); }

	CooIterator operator+(int n) const
		{ return CooIterator(parent, i+n); }
	bool operator==(CooIterator const &rhs) const
		{ return i == rhs.i; }
	bool operator!=(const CooIterator& rhs) const
		{return !this->operator==(rhs); }


	int offset() const { return i; }
	IterIndexT &index(int k)
		{ return parent->index(k,i); }
	void set_index(indices_type const &idx)
		{ parent->set_index(i, idx); }
	IterValT &val() { return parent->val(i); }
};
// -----------------------------------------------------

template<class ArrayT>
std::ostream &_ostream_out_array(std::ostream &os, ArrayT const &A);

template<class ArrayT>
std::ostream &_ostream_out_array(std::ostream &os, ArrayT const &A)
{
	os << "VectorCooArray<";
	stream(os, &A.shape[0], A.shape.size());
	os << ">(";
	for (auto ii(A.begin()); ii != A.end(); ++ii) {
		os << "(";
		auto idx(ii.index());
		for (int k=0; k<A.rank; ++k) os << idx[k] << " ";
		os << ": " << ii.val() << ")";
	}
	os << ")";
	return os;
}

/** @} */

}	// Namespace



// ---------------------------------------------------------------------------


#endif	// Guard
