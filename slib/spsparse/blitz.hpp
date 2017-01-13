#ifndef SPSPARSE_BLITZ_HPP
#define SPSPARSE_BLITZ_HPP

#include <ibmisc/ibmisc.hpp>
#include <spsparse/spsparse.hpp>
#include <ibmisc/blitz.hpp>

namespace spsparse {

#if 0        // see ibmisc/blitz.hpp
template<class AccumulatorT, class ValT, int RANK>
extern void spcopy(AccumulatorT &ret, blitz::Array<ValT, RANK> const &A, bool set_shape=true);

template<class AccumulatorT, class ValT, int RANK>
void spcopy(AccumulatorT &ret, blitz::Array<ValT, RANK> const &A, bool set_shape)
{
    if (set_shape) ret.set_shape(A.shape);
    std::array<int,AccumulatorT::rank> idx;
    for (auto ii=A.begin(); ii != A.end(); ++ii) {
        ret.add(ii.index(), ii.val());
    }
}
#endif

// =========================================================

// Accumulator for use with Spsparse
/** @brief For output/conversion to dense arrays.

Outputs to a blitz::Array.  Note that blitz:Array is quite flexible,
and can be used to access almost any existing fixed-size data
structure.

Usage Example:
@code
VectorCooArray<int,double,2> A;
blitz::Array<double,2> B(to_tiny(A.shape));
BlitzAccum<decltype(A)> Baccum(B);
spcopy(Baccum, A);
@endcode

*/
template<class ValueT, int RANK>
struct BlitzAccum
{
    static const int rank = RANK;
    typedef int index_type;
    typedef ValueT value_type;

private:
    blitz::Array<value_type,rank> * const result;
    DuplicatePolicy duplicate_policy;
    bool shape_is_set;
    double fill_value;

public:
    BlitzAccum(
        blitz::Array<value_type,rank> * const _dense,
        bool reset_shape,    // Should we re-allocate on set_shape() call?
        ValueT _fill_value,
        DuplicatePolicy _duplicate_policy)
    : result(_dense), shape_is_set(!reset_shape),
        fill_value(_fill_value), duplicate_policy(_duplicate_policy)
    {
    }

    void set_shape(std::array<long, RANK> const &_shape);

    inline void add(std::array<index_type,rank> const &index, value_type const &val);
};

template<class ValueT, int RANK>
void BlitzAccum<ValueT,RANK>::set_shape(std::array<long, RANK> const &_shape)
{
    if (!shape_is_set) {
        blitz::TinyVector<int,rank> shape_t;
        for (int i=0; i<RANK; ++i) shape_t[i] = _shape[i];
        result->reference(blitz::Array<value_type,rank>(shape_t));
        *result = fill_value;
        shape_is_set = true;
    }
}


template<class ValueT, int RANK>
inline void BlitzAccum<ValueT,RANK>::add(std::array<index_type,rank> const &index, value_type const &val)
{
    // Check bounds...
    blitz::TinyVector<int, rank> bidx;
    for (unsigned int i=0; i<rank; ++i) {
        auto &ix(index[i]);
        if (ix < result->lbound(i) || ix > result->ubound(i)) (*ibmisc::ibmisc_error)(-1,
            "Index %d out of bounds: %d vs [%d, %d]",
            i, ix, result->lbound(i), result->ubound(i));
        bidx[i] = index[i];
    }
    value_type &oval((*result)(bidx));

    switch(duplicate_policy) {
        case DuplicatePolicy::LEAVE_ALONE :
            if (!std::isnan(oval)) oval = val;
        break;
        case DuplicatePolicy::ADD :
            oval += val;
        break;
        case DuplicatePolicy::REPLACE :
            oval = val;
        break;
        case DuplicatePolicy::REPLACE_THEN_ADD :
            if (std::isnan(oval)) oval = val;
            else oval += val;
        break;
    }
}




template<class ValueT, int RANK>
inline BlitzAccum<ValueT, RANK> blitz_accum_new(
    blitz::Array<ValueT, RANK> * const _dense,
    ValueT _fill_value=0,
    DuplicatePolicy _duplicate_policy = DuplicatePolicy::ADD)
{ return BlitzAccum<ValueT,RANK>(_dense, true, _fill_value, _duplicate_policy); }

template<class ValueT, int RANK>
inline BlitzAccum<ValueT, RANK> blitz_accum_existing(
    blitz::Array<ValueT, RANK> * const _dense,
    DuplicatePolicy _duplicate_policy = DuplicatePolicy::ADD)
{ return BlitzAccum<ValueT,RANK>(_dense, false, 0, _duplicate_policy); }

// ----------------------------------------------------------

template<class AccumulatorT, class TypeT, int RANK>
extern void spcopy(AccumulatorT &ret,
    blitz::Array<TypeT, RANK> const &arr, bool set_shape=true);

template<class AccumulatorT, class TypeT, int RANK>
void spcopy(AccumulatorT &ret,
    blitz::Array<TypeT, RANK> const &arr, bool set_shape)
{
    if (set_shape) {
        std::array<long,RANK> shape;
        for (int i=0; i<RANK; ++i) shape[i] = arr.extent(i);
        ret.set_shape(shape);
    }

    typedef typename AccumulatorT::index_type IndexT;
    for (auto ii=arr.begin(); ii != arr.end(); ++ii) {
        if (*ii != 0) {
            auto index(ibmisc::to_array<IndexT,int,RANK>(ii.position()));
            ret.add(index, *ii);
        }
    }
}
// ----------------------------------------------------------
/** Convert spsparse-object to blitz, as long as it has an associated spcopy() function. */
template<class SourceT>
extern blitz::Array<typename SourceT::value_type, SourceT::rank>
to_blitz(SourceT const &M);

template<class SourceT>
extern blitz::Array<typename SourceT::value_type, SourceT::rank>
to_blitz(SourceT const &M)
{
    blitz::Array<typename SourceT::value_type, SourceT::rank> ret;
    auto accum1(blitz_accum_new(&ret));
    spcopy(accum1, M, true);
    return ret;
}



}    // namespace

#endif
