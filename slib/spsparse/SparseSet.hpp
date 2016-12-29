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

#include <unordered_map>
#include <algorithm>
#include <spsparse/spsparse.hpp>
#include <ibmisc/array.hpp>
#include <spsparse/blitz.hpp>

namespace spsparse {

/** Translates between a sparse set (say, the set of indices used in a
SparseVector) and a dense set numbered [0...n) */
template<class SparseT, class DenseT>
class SparseSet {
    SparseT _sparse_extent;
    std::unordered_map<SparseT, DenseT> _s2d;
    std::vector<SparseT> _d2s;

public:
    SparseSet() : _sparse_extent(-1) {}

    void clear() {
        _sparse_extent = -1;
        _s2d.clear();
        _d2s.clear();
    }

    bool in_sparse(SparseT const &sparse_ix) const
        { return _s2d.find(sparse_ix) != _s2d.end(); }

    bool in_dense(DenseT dense_ix) const
        { return (dense_ix >= 0 && dense_ix < dense_extent()); }

    SparseT sparse_extent() const
        { return _sparse_extent; }

    void set_sparse_extent(SparseT extent)
    {
        std::stringstream buf;
        if (_sparse_extent != -1 && _sparse_extent != extent) {
            buf << "Cannot change sparse_extent from " << _sparse_extent << " to " << extent;
            (*ibmisc::ibmisc_error)(-1, "%s", buf.str().c_str());
        }
        _sparse_extent = extent;
    }

    DenseT dense_extent() const
        { return _d2s.size(); }

    void add(SparseT sparse_index)
    {
        if (_s2d.find(sparse_index) == _s2d.end()) {
            _s2d.insert(std::pair<SparseT,DenseT>(sparse_index, dense_extent()));
            _d2s.push_back(sparse_index);
        }
    }

    template<class IterT>
    void add(IterT sparsei, IterT const &sparse_end)
    {
        for (; sparsei != sparse_end; ++sparsei) add(*sparsei);
    }


    template<class IterT>
    void add_sorted(IterT sparsei, IterT const &sparse_end)
    {
        std::vector<SparseT> sorted;
        for (; sparsei != sparse_end; ++sparsei) sorted.push_back(*sparsei);
        std::sort(sorted.begin(), sorted.end());

        add(sorted.begin(), sorted.end());
    }

    DenseT to_dense(SparseT const &sval) const
        { return _s2d.at(sval); }

    SparseT to_sparse(DenseT const &dval) const
    {
        if (dval < 0 || dval >= _d2s.size()) {
            (*ibmisc::ibmisc_error)(-1,
                "Value %ld is out of range (0, %ld)", dval, _d2s.size());
        }
        return _d2s[dval];
    }
};

// -----------------------------------------------------------
enum class SparseTransform {ID, TO_DENSE, TO_SPARSE};

template<class AccumulatorT,class SparseSetT>
class SparseTransformAccum
{
    struct Data {
        SparseSetT const * const sparse_set;
        SparseTransform transform;

        Data(SparseSetT *_sparse_set, SparseTransform _transform) :
            sparse_set(_sparse_set), transform(_transform) {}
    };

public:
    SPSPARSE_LOCAL_TYPES(AccumulatorT);

private:
    AccumulatorT * const sub;
    std::vector<Data> data;

public:
    /** @param _transform
        SparseTransform::TO_DENSE
        SparseTransform::TO_SPARSE
    @param sparse_sets
        0 if that dimension is not to be transformed
    */
    SparseTransformAccum(
        AccumulatorT *_sub,
        SparseTransform _transform,
        std::array<SparseSetT *, (size_t)rank> const sparse_sets);

    /** Merge shape we're given with shape from transforms. */
    void set_shape(std::array<long, rank> shape)
    {
        for (int i=0; i<rank; ++i) {
            // Use the shape we were given, if no transform for this dimension
            if (!data[i].sparse_set) continue;

            switch(data[i].transform) {
                case SparseTransform::TO_DENSE:
                    shape[i] = data[i].sparse_set->dense_extent();
                    break;
                case SparseTransform::TO_SPARSE:
                    shape[i] = data[i].sparse_set->sparse_extent();
                    break;
            }
        }
        sub->set_shape(shape);
    }

    void add(std::array<index_type,rank> index, val_type const &val)
    {
        for (int i=0; i<rank; ++i) {
            // Use the index we were given, if no transform for this dimension
            if (!data[i].sparse_set) continue;

            switch(data[i].transform) {
                case SparseTransform::TO_DENSE:
                    index[i] = data[i].sparse_set->to_dense(index[i]);
                    break;
                case SparseTransform::TO_SPARSE:
                    index[i] = data[i].sparse_set->to_sparse(index[i]);
                    break;
            }
        }
        sub->add(index, val);
    }
};

template<class AccumulatorT,class SparseSetT>
SparseTransformAccum<AccumulatorT,SparseSetT>::SparseTransformAccum(
    AccumulatorT *_sub,
    SparseTransform transform,
    std::array<SparseSetT *, (size_t)rank> const sparse_sets)    // 0 if we don't want to densify/sparsify that dimension
    : sub(_sub)
{
    for (int i=0; i<rank; ++i) {
        data.push_back(Data(sparse_sets[i],
            sparse_sets[i] ? transform : SparseTransform::ID));
    }
}

#define SparseTransformAccumT SparseTransformAccum<AccumulatorT,SparseSetT>
template<class AccumulatorT,class SparseSetT>
inline SparseTransformAccumT sparse_transform_accum(
        AccumulatorT *sub,
        SparseTransform transform,
        std::array<SparseSetT *, AccumulatorT::rank> const &sparse_sets)
{
    return SparseTransformAccumT(sub, transform, sparse_sets);
}
#undef SparseTransformAccumT



template<class AccumulatorT, class SrcT, class SparseT, class DenseT>
void sparse_copy(AccumulatorT &ret, SrcT const &src,
    SparseTransform direction,
    std::array<spsparse::SparseSet<SparseT, DenseT> *, AccumulatorT::rank> dims,
    bool set_shape=true);

template<class AccumulatorT, class SrcT, class SparseT, class DenseT>
void sparse_copy(AccumulatorT &ret, SrcT const &src,
    SparseTransform direction,
    std::array<spsparse::SparseSet<SparseT, DenseT> *, AccumulatorT::rank> dims,
    bool set_shape)
{
    auto accum(sparse_transform_accum(&ret, direction, dims));
    spcopy(accum, src, set_shape);
}




template<class AccumulatorT, class SrcT, class SparseT, class DenseT>
void sparse_copy(AccumulatorT &ret, SrcT const &src,
    SparseTransform direction,
    std::array<spsparse::SparseSet<SparseT, DenseT> const *, AccumulatorT::rank> dims,
    bool set_shape=true);

template<class AccumulatorT, class SrcT, class SparseT, class DenseT>
void sparse_copy(AccumulatorT &ret, SrcT const &src,
    SparseTransform direction,
    std::array<spsparse::SparseSet<SparseT, DenseT> const *, AccumulatorT::rank> dims,
    bool set_shape)
{
    auto accum(sparse_transform_accum(&ret, direction, dims));
    spcopy(accum, src, set_shape);
}


// -------------------------------------------------------
/** A data sructure that acts like a normal SpSparse Accumulator.  BUT:
  1) When stuff is added to it, it also updates corresponding dimension maps (SparseSet).
  2) An Eigen::SparseMatrix may be extracted from it when construction is complete. */
template <class AccumulatorT, class DenseIndexT>
class MappedArray
{
public:
    typedef typename AccumulatorT::index_type SparseIndexT;
    typedef typename AccumulatorT::val_type ValT;
    typedef SparseSet<SparseIndexT, DenseIndexT> SparseSetT;
    static const int rank = AccumulatorT::rank;

    /** The underlying constructed matrix.  Access directly if needed. */
    AccumulatorT M;
    std::array<SparseSetT *, rank> const dims;

    /** @param dims Dimension maps for each dimension.  If one is
    preparing to multiply matrices, then each dimension map will be
    shared by at least two MappedArray object. */
    MappedArray(std::array<SparseSetT *, rank> const &_dims) : dims(_dims) {}

    /** Sets the shape of the underlying matrix (and dimension maps) */
    void set_shape(std::array<long, rank> const &shape)
    {
        M.set_shape(shape);
        for (int k=0; k<rank; ++k) dims[k]->set_sparse_extent(shape[k]);
    }

    /** Adds an item. */
    void add(std::array<SparseIndexT, rank> const &index, ValT const val)
    {
        M.add(index, val);
        for (int k=0; k<rank; ++k) dims[k]->add(index[k]);
    }

};

// -----------------------------------------------------------




}   // Namespace
