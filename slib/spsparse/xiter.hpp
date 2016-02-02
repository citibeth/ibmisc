#ifndef XITER_HPP
#define XITER_HPP

#include <type_traits>

namespace spsparse {

/** @defgroup xiter xiter.hpp
@brief Extended iterator types used by SpSparse

The fundamental type of iterator used for sparse-sparse matrix multiplation is the join iterator.  It takes two (or more) SORTED sequences of elements, and returns values only when all sequences match value.  For example:
@code
    join_iterator([0, 3, 4, 8], [1, 4, 5, 6, 7, 8, 10])
@endcode
prodces the values:
@code
    [1,4,8]
@endcode
See spsparse::Join2Xiter and spsparse::Join3Xiter for variants that
join two and three sub-iterators.


Two types of extensions to the usual STL iterator interface were needed for SpSparse:

1. In Python, when iterating through a dict (a set of key/value
   pairs), the "default" iterator iterates through keys.  Similarly,
   when iterating through a sparse matrix, the "default" iterator
   should return the "key" (std::array<int,RANK> index) of each item
   in the matrix.  The method val() is used to return the value of an
   element.  For example:
@code
VectorCooMatrix<int, double> A;
for (auto ii = A.begin(); ii != A.end(); ++ii) {
    std::array<int,2> index(*ii);
    printf("A[%d,%d] = %g", index[0], index[1], ii.val();
}
@endcode

   We call an iterator conforming to this property a ValIterator.
   Note that the iterators on VectorCooArray also add other ways to access
   the index, but these are not part of the ValIterator interface:
@code
VectorCooMatrix<int, double> A;
for (auto ii = A.begin(); ii != A.end(); ++ii)
    printf("A[%d,%d] = %g", ii.index(0), ii.index(1), ii.val();
@endcode
2. Iterators that internally know their own end (Xiter).  See
   spsparse::STLXiter for a way to wrap any STL forward iterator to
   conform to this interface. See also spsparse::ValSTLXiter


@{
*/


// -----------------------------------------------
// -----------------------------------------------
/** @brief Convert standard STL iterator into our XIter.

Code Example:
@code
	std::vector<int> vec;
	for (auto ii(make_xiter(vec.begin(), vec.end());
		!ii.eof(); ++ii)
	{ printf("vec[%d] = %d\n", ii.offset(), *ii); }
@endcode

*/
template<class STLIter>
class STLXiter {
public:
	/** @brief Beginning of iteration. */
	STLIter const begin;
	/** @brief Current position of iteration (access as needed). */
	STLIter ii;
	/** @brief End of iteration */
	STLIter const end;

	typedef typename STLIter::value_type value_type;

	STLXiter(STLIter const &_begin, STLIter const &_end) :
		begin(_begin), ii(_begin), end(_end) {}

	/** @brief True if iterator has no more values (cannot be dereferenced). */
	bool eof() { return (ii == end); }

	/** @brief Number of times the iterator has been incremented
	since the beginning. */
	size_t offset() { return ii - begin; }

	/** @brief Move the iterator forward one step. */
	void operator++() { ++ii; }

	/** @brief Pass-through, dereference the iterator. */
	auto operator*() -> decltype(*ii) { return *ii; }
};

/** @brief Convert standard STL iterator (with extra .val() method)
into our XIter.

Code Example:
@code
	VectorCooVector<int, double> vec;
	for (auto ii(make_val_xiter(vec.begin(), vec.end());
		!ii.eof(); ++ii)
	{ printf("vec[%d] = %d\n", ii.index(0), ii.val()); }
@endcode

@see spsparse::STLXiter */
template<class ValSTLIter>
class ValSTLXiter : public STLXiter<ValSTLIter>
{
public:

	ValSTLXiter(ValSTLIter const &_begin, ValSTLIter const &_end) :
		STLXiter<ValSTLIter>(_begin, _end) {}

	/** @brief Pass-through, call val() on the iterator. */
	auto val() -> decltype(STLXiter<ValSTLIter>::ii.val())
		{ return this->ii.val(); }
};
// --------------------------------------------------------

// --------------------------------------------------------

// http://stackoverflow.com/questions/984394/why-not-infer-template-parameter-from-constructor

template<class STLIter>
STLXiter<STLIter> make_xiter(
	STLIter **_begin, STLIter &&_end)
{ return STLXiter<STLIter>(std::move(_begin), std::move(_end)); }

template<class ValSTLIter>
ValSTLXiter<ValSTLIter> make_val_xiter(
	ValSTLIter &&_begin, ValSTLIter &&_end)
{ return ValSTLXiter<ValSTLIter>(std::move(_begin), std::move(_end)); }


// -------------------------------------------------------





// -------------------------------------------------------
/** @brief Joins three (sorted, non-repeating) Xiters.  See Join2Xiter.
@see Join2Xiter
*/
template<class Xiter1T, class Xiter2T, class Xiter3T>
class Join3Xiter
{
	typename std::remove_const<typename Xiter1T::value_type>::type next_match;
	bool _eof;

public:
	// Allow user to access underlying Xiters, to get useful stuff out of them.
	Xiter1T i1;
	Xiter2T i2;
	Xiter3T i3;
	int total_in_use;

	bool eof() { return _eof; }

	Join3Xiter(Xiter1T &&_i1, Xiter2T &&_i2, Xiter3T &&_i3) :
		next_match(0),
		i1(std::move(_i1)),
		i2(std::move(_i2)),
		i3(std::move(_i3))
	{
		_eof = i1.eof();
		if (_eof) return;
		next_match = *i1;
		next_noincr();
	}

private:
	void next_noincr()
	{
#define JOIN_RANK 3
#include "next_noincr_body.hpp"
#undef JOIN_RANK
	}

public:
	void operator++()
	{
		++i1;
		++i2;
		++i3;
		// We now don't know, i1 or i2 might now be EOF
		next_noincr();
	}

};

template<class Xiter1T, class Xiter2T, class Xiter3T>
class Join3Xiter<Xiter1T, Xiter2T, Xiter3T> join3_xiter(Xiter1T &&i1, Xiter2T &&i2, Xiter3T &&i3)
	{ return Join3Xiter<Xiter1T, Xiter2T, Xiter3T>(std::move(i1), std::move(i2), std::move(i3)); }

// ========================================================
/** @brief Joins two (sorted, non-repeating) Xiters, producing output when the two match.

Code Example:
@code
std::vector<int> v1 = {0, 3, 4, 8};
std::vector<int> v2 = {1, 4, 5, 6, 7, 8, 10};
for (auto ii(join2_xiter(
	make_xiter(v1.begin(), v1.end()),
	make_xiter(v2.begin(), v2.end())));
	!ii.eof(); ++ii)
{
	printf("Matching value is %d\n", *v1.i1);
}
@endcode

@note This Xiter does NOT directly provide a dereference operator.  To
fish values out of your iterator, access the sub-iterators directly.
For example:

@code
VectorCooVector<int, double> v1, v2;
for (auto ii(join2_xiter(
	make_xiter(v1.begin(), v1.end()),
	make_xiter(v2.begin(), v2.end())));
	!ii.eof(); ++ii)
{
	// NOTE: *ii.v1 == *ii.v2 here!
	printf("v1[%d] = %g, v2[%d] = %g\n",
		*ii.v1, ii.v1.val(), *ii.v2, ii.v2.val());
}
@endcode


@see spsparse::join2_xiter
*/
template<class Xiter1T, class Xiter2T>
class Join2Xiter
{
	typename std::remove_const<typename Xiter1T::value_type>::type next_match;
	bool _eof;

public:
	// Allow user to access underlying Xiters, to get useful stuff out of them.
	Xiter1T i1;
	Xiter2T i2;
	int total_in_use;

	bool eof()
		{ return _eof; }

	Join2Xiter(Xiter1T &&_i1, Xiter2T &&_i2) :
		i1(std::move(_i1)),
		i2(std::move(_i2))
	{
		_eof = i1.eof();
		if (_eof) return;
		next_match = *i1;
		next_noincr();
	}

private:
	void next_noincr()
	{
#define JOIN_RANK 2
#include "next_noincr_body.hpp"
#undef JOIN_RANK
	}

public:
	void operator++()
	{
		++i1;
		++i2;
		// We now don't know, i1 or i2 might now be EOF
		next_noincr();
	}

};

template<class Xiter1T, class Xiter2T>
class Join2Xiter<Xiter1T, Xiter2T> join2_xiter(Xiter1T &&i1, Xiter2T &&i2)
	{ return Join2Xiter<Xiter1T, Xiter2T>(std::move(i1), std::move(i2)); }

// ========================================================

/** @} */

}	// Namespace

#endif
