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

#ifndef SPSPARSE_SPSPARSE_H
#define SPSPARSE_SPSPARSE_H

#include <exception>
#include <cmath>
#include <ibmisc/iter.hpp>

namespace spsparse {

/** @defgroup spsparse spsparse.hpp
@brief Basic stuff common to all of SpSparse library.

@{
*/


/** @brief What to do in algorithms when duplicate entries are encountered.

- ADD (default): Sum them together.
- LEAVE_ALONE: Use the first value encountered.
- REPLACE: Use the last value encountered.
- REPLACE_THEN_ADD: (for dense destinations only) Add if it's not NaN, otherwise replace

@see spsparse::consolidate() */
enum class DuplicatePolicy {
    LEAVE_ALONE, ADD, REPLACE, REPLACE_THEN_ADD};


/** @brief Promote relevant template parameters.

Used to propagate relevant template parameters throughout the
different classes that need them.  Provides: rank, index_type, val_type.

@note val_type is DIFFERENT from the standard STL val_type.  Standard STL val_type is the same as our index_type.

Code Example
@code
template<class VectorCooArrayT>
class MyClass {
public:
    SPSPARSE_LOCAL_TYPES(VectorCooArrayT);

};
@endcode
*/
template<class ArrayOrIterT>
struct AccumTraits {
    static const int rank = ArrayOrIterT::rank;
    typedef typename ArrayOrIterT::index_type index_type;
    typedef typename ArrayOrIterT::val_type val_type;
};



// Deprecated...
#define SPSPARSE_LOCAL_TYPES(ArrayOrIterT) \
    static const int rank = ArrayOrIterT::rank; \
    typedef typename ArrayOrIterT::index_type index_type; \
    typedef typename ArrayOrIterT::val_type val_type;

// -------------------------------------------------------------

// The expression "std::isnan(n) || (n == 0)" for different data types.
// Use template specilization here...
/** @brief Internal helper function.

Used to tell if a value is "none" (i.e. 0 or NaN) and should therefore
be eliminated under zero_nan.  This will need to be specialized for
value types that are not a single IEEE floating point (eg:
std::complex<double>).

@see spsparse::consolidate() */
template<class NumT>
inline bool isnone(NumT const n, bool const zero_nan=false)
{
    if (zero_nan) {
        return std::isnan(n) || (n == 0);
    } else {
        return (n == 0);
    }
}



// -----------------------------------------------------
/** @brief Select out just one dimension of the index on iteration.

Wraps EigenSparseMatrixIterator or TupleList::iterator, producing a
single dimension's index with operator*()

Code Example
@code
VectorCooMatrix<int, double> const A;
typedef DimIndexIter<decltype(A)::const_iterator> DIType;
for (DIType ii(1, A.begin()); ii != DIType(1, A.end()); ++ii)
    printf("Element with column %d and value %d\n", *ii, ii.val());
@endcode

@see spsparse::VectorCooMatrix::dim_iter(), spsparse::VectorCooMatrix::dim_begin(), spsparse::VectorCooMatrix::dim_end()
*/
template<class IterT>
class DimIndexIter : public ibmisc::forward_iterator<typename IterT::index_type, DimIndexIter<IterT>>
{
public:
    // Type traits for spsparse iterators
    static int const rank = IterT::rank;
    typedef typename IterT::index_type index_type;
    typedef typename IterT::val_type val_type;
private:
    IterT sub;
    const int dim;

public:
    DimIndexIter(IterT &&ii, int _dim) : sub(std::move(ii)), dim(_dim) {}

    val_type operator*()
        { return sub->index(dim); }

    DimIndexIter &operator++()
        { ++sub; return *this; }
    bool operator==(const DimIndexIter& rhs) const
        {return sub == rhs.sub;}
};

template<class IterT>
DimIndexIter<IterT> dim_index_iter(IterT &&sub, int dim)
    { return DimIndexIter<IterT>(std::move(sub), dim); }
// ------------------------------------------------------------------------




} // Namespace spsparse

/** @} */


#endif // Guard
