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
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/tuple.hpp>

namespace spsparse {
namespace accum {

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
template<class AccumT>
class Filter : public AccumTraits<AccumT>
{
    typedef AccumTraits<AccumT> super;
public:
    AccumT sub;

    Filter(AccumT &&_sub)
        : sub(std::move(_sub)) {}

    typename super::base_array_type &base()
        { return sub.base(); }

    void set_shape(std::array<long, super::rank> const &_shape)
        { sub.set_shape(_shape); }

    void add(
        std::array<typename super::index_type, super::rank> const &index,
        typename super::val_type const &val)
        { sub.add(index, val); }
};
// -----------------------------------------------------------
/* Like the Unix "tee", sends output to more than one accumulator
destinations.  For example:

    spcopy(refs(accum_dest1, accum_dest2), accum_src);

The number of destinations is arbitrary.
*/
template<class AccumT0, class... AccumTs>
class Tee : public AccumTraits<AccumT0>
{
    typedef AccumTraits<AccumT0> super;
    boost::fusion::tuple<AccumT0, AccumTs...> subs;
public:

    Tee(AccumT0 &&_sub0, AccumTs&&... _subs)
        : subs(std::move(_sub0), std::move(_subs)...) {}

    typename super::base_array_type &base()
        { return boost::fusion::get<0>(subs).base(); }


    // ----------------------------------------------
private:
    // Helper class for set_shape().  See
    //    https://theboostcpplibraries.com/boost.fusion
    struct _set_shape
    {
        std::array<long, super::rank> const &shape;

        _set_shape(std::array<long, super::rank> const &_shape)
        : shape(_shape) {}

        template<class typei>
        void operator()(typei &sub) const
            { sub.set_shape(shape); }
    };
public:
    void set_shape(std::array<long, super::rank> const &shape)
    {
        boost::fusion::for_each(subs, _set_shape(shape));
    }
    // ----------------------------------------------
private:
    // Helper class for add().  See
    //    https://theboostcpplibraries.com/boost.fusion
    struct _add
    {

        std::array<typename super::index_type, super::rank> const &index;
        typename super::val_type const &val;

        _add(std::array<typename super::index_type, super::rank> const &_index,
            typename super::val_type const &_val)
        : index(_index), val(_val) {}

        template<class typei>
        void operator()(typei &sub) const
            { sub.add(index, val); }
    };
public:
    void add(
        std::array<typename super::index_type, super::rank> const &index,
        typename super::val_type const &val)
    {
        boost::fusion::for_each(subs, _add(index, val));
    }
};

template<class AccumT0, class ...AccumTs>
Tee<AccumT0, AccumTs...> tee(AccumT0 &&sub0, AccumTs&&... subs)
    { return Tee<AccumT0, AccumTs...>(std::move(sub0), std::move(subs)...); }
// -----------------------------------------------------------
/** Accumulator filters generally own the sub-accumulator upon which
they are stacked; this makes it easy to create them on-the-fly in a
spcopy() call.  However, these accumulators disappear once the
spcopy() call is over, providing no way to accumulate into something
permanent.  This shim references an already-created filter, allowing
other filtes to wrap it on-the-fly.  For example:

    // Simple copy is OK
    TupleList<...> accum_dest;
    spcopy(accum_dest, accum_src);

    // Wrapping a raw destination is not OK.
    spcopy(transposer(accum_dest), accum_src);

    // Wrapping a raw destination is not OK; the transpose() filter needs
    // to own its sub-filter.
    spcopy(transposer(accum_dest), accum_src);

    // This make it work:
    spcopy(transposer(ref(accum_dest)), accum_src);
*/
template<class AccumT>
class Ref : public AccumTraits<AccumT>
{
    typedef AccumTraits<AccumT> super;
public:
    AccumT &sub;

    Ref(AccumT &_sub)
        : sub(_sub) {}

    typename super::base_array_type &base()
        { return sub.base(); }

    void set_shape(std::array<long, super::rank> const &_shape)
        { sub.set_shape(_shape); }

    void add(
        std::array<typename super::index_type, super::rank> const &index,
        typename super::val_type const &val)
        { sub.add(index, val); }
};

template<class AccumT>
Ref<AccumT> ref(AccumT &sub)
    { return Ref<AccumT>(sub); }
// -----------------------------------------------------------
/** @brief For transpose or project.

Permutes/selects dimensions, storing in a sub-accumulator.  The
sub-accumulator does NOT have to have the same rank; this may
therefore be used to transpose or remove dimensions.

Usage Example:
@code
VectorCooArray<int,double,2> A,B;
PermuteAccum<A::rank, decltype(B)> p(B, {1,0});
spcopy(p, A);
@endcode

*/
template<class AccumT,size_t IN_RANK>
class Permute : public Filter<AccumT>
{
    typedef Filter<AccumT> super;
public:
    // Override rank stuff from FilterAccum
    static const size_t rank = IN_RANK;
    static const size_t out_rank = AccumT::rank;

private:
    std::array<int,out_rank> perm;
    std::array<int, out_rank> out_idx;

public:
    Permute(std::array<int,out_rank> const _perm, AccumT &&_sub)
        : super(std::move(_sub)), perm(_perm) {}

    void set_shape(std::array<long, rank> const &_shape)
    {
        std::array<long, out_rank> oshape;
        for (int i=0; i<out_rank; ++i) oshape[i] = _shape[perm[i]];
        super::sub.set_shape(oshape);
    }

    void add(
        std::array<typename super::index_type,rank> const &index,
        typename super::val_type const &val)
    {
        for (int i=0; i<out_rank; ++i) out_idx[i] = index[perm[i]];
        super::sub.add(out_idx, val);
    }
};

// Helper used to specify template arguments
template<size_t X>
struct in_rank {
    static const size_t x=X;
};

template<class AccumT,size_t IN_RANK>
inline Permute<AccumT,IN_RANK>
permute(
    in_rank<IN_RANK> const in_rank_dummy,
    std::array<int, AccumT::rank> const &_perm,
    AccumT &&_sub)
    { return Permute<AccumT,IN_RANK>(_perm, std::move(_sub)); }

// -----------------------------------------------------------
template<class AccumT>
class Transpose : public Permute<AccumT,2>
{
public:
    Transpose(AccumT &&_sub) :
        Permute<AccumT,2>({1,0}, std::move(_sub)) {}
};


template<class AccumT>
Transpose<AccumT>
transpose(AccumT &&_sub)
    { return Transpose<AccumT>(std::move(_sub)); }

// -----------------------------------------------------------
template<class AccumT, class TransformFn>
class TransformAccum : public Filter<AccumT>
{
    typedef Filter<AccumT> super;
    TransformFn transform_fn;
public:

    TransformAccum(TransformFn &&_transform_fn, AccumT &&_sub)
        : super(std::move(_sub)), transform_fn(std::move(_transform_fn)) {}
    void add(
        std::array<typename super::index_type,super::rank> const &index,
        typename super::val_type const &val)
    {
        super::sub.add(index, transform_fn(val));
    }

};

// -----------------------------------------------------------
struct InvertFn{
    char invert;
    InvertFn(char _invert) : invert(_invert) {}
    double operator()(double x)
        { return invert=='-' ? 1./x : x; }
};
template<class AccumT>
TransformAccum<AccumT,InvertFn> invert(char invert, AccumT &&sub)
    { return TransformAccum<AccumT,InvertFn>(InvertFn(invert), std::move(sub)); }
// -------------------------------------------------------




/** @} */


}}   // Namespace

#endif // Guard
