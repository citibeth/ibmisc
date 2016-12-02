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

#ifndef SPSPARSE_ACCUM_HPP
#define SPSPARSE_ACCUM_HPP

#include <cstddef>
#include <array>
#include <vector>
#include <ibmisc/blitz.hpp>
#include <spsparse/spsparse.hpp>

namespace spsparse {

/** @defgroup accum accum.hpp
@brief Accumulators to use with algorithms.

An accumulator is used by an algorithm to accept output via an add()
method.  This allows each algorithm to support a variety of output
strategies: in-place, out-of-place, conversion to dense arrays, sum
for inner product, etc.

TODO: Figure out an implement good semantics for set_shape().  When is
it called / not called?

@{
*/

// -----------------------------------------------------------
/** @brief For in-place operations.

Overwrites a VectorCooArray::iterator.  Good for
in-place operations that we KNOW won't change the size of the
VectorCooArray (eg: use with spsparse::transpose()).

@warning Does not bounds-check for the END of the iterator.

Usage Example:
@code
VectorCooArray<int,double,2> A;
OverWriteAccum<decltype(A)::iterator> overwrite(A.begin());
transpose(overwrite, A, {1,0});
@endcode

*/
template<class IterT>
class OverwriteAccum
{
    SPSPARSE_LOCAL_TYPES(IterT);

    IterT ii;
public:
    OverwriteAccum(IterT &&_ii) : ii(std::move(_ii)) {}

    void add(std::array<index_type,rank> const &index, typename IterT::val_type const &val) {
        ii.set_index(index);
        ii.val() = val;
        ++ii;
    }
};
// -----------------------------------------------------------
/** @brief For transpose or project.

Permutes/selects dimensions, storing in a sub-accumulator.  The
sub-accumulator does NOT have to have the same rank; this may
therefore be used to transpose or remove dimensions.

Usage Example:
@code
VectorCooArray<int,double,2> A,B;
PermuteAccum<A::rank, decltype(B)> p(B, {1,0});
copy(p, A);
@endcode

*/
template<class AccumulatorT,size_t IN_RANK>
class PermuteAccum
{
public:
    static const size_t rank = IN_RANK;
    static const size_t out_rank = AccumulatorT::rank;

    typedef typename AccumulatorT::index_type index_type;
    typedef typename AccumulatorT::val_type val_type;

private:
    AccumulatorT * const sub;
    std::array<int,out_rank> perm;
    std::array<int, out_rank> out_idx;

public:
    PermuteAccum(AccumulatorT *_sub, std::array<int,out_rank> const &_perm)
        : sub(_sub), perm(_perm) {}

    void set_shape(std::array<long, rank> const &_shape)
    {
        std::array<long, out_rank> oshape;
        for (int i=0; i<out_rank; ++i) oshape[i] = _shape[perm[i]];
        sub->set_shape(oshape);
    }

    void add(std::array<index_type,rank> const &index, val_type const &val) {
        for (int i=0; i<out_rank; ++i) out_idx[i] = index[perm[i]];
        sub->add(out_idx, val);
    }
};
// -----------------------------------------------------------
template<size_t X>
struct in_rank {
    static const size_t x=X;
};

template<class AccumulatorT,size_t IN_RANK>
inline PermuteAccum<AccumulatorT,IN_RANK>
permute_accum(
    AccumulatorT *_sub,
    in_rank<IN_RANK> const in_rank_dummy,
    std::array<int, AccumulatorT::rank> const &_perm)
    { return PermuteAccum<AccumulatorT,IN_RANK>(_sub, _perm); }

// -----------------------------------------------------------

// -----------------------------------------------------------


// -------------------------------------------------------
/** @brief For inner product.

Accumulates into a single scalar, ignoring index information.  Used to
implement inner products, or just sum over an array.

Usage Example:
@code
VectorCooArray<int,double,2> A;
ScalarAccumulator<decltype(A)> s;
copy(s, A);
printf("Sum = %g\n", s.val);
@endcode
*/
template<class VectorCooArrayT>
struct ScalarAccumulator {
    SPSPARSE_LOCAL_TYPES(VectorCooArrayT);

    val_type val;
    ScalarAccumulator() : val(0) {}

    void add(const std::array<index_type, rank> &index, val_type const &_val)
        { this->val += _val; }
};




/** @} */


}   // Namespace

#endif // Guard
