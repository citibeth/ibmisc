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

private:
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

public:

    template<class IterT>
    void add_sorted(IterT sparsei, IterT const &sparse_end)
    {
        std::vector<SparseT> sorted;
        for (; sparsei != sparse_end; ++sparsei) sorted.push_back(*sparsei);
        std::sort(sorted.begin(), sorted.end());

        add(sorted.begin(), sorted.end());
    }

    DenseT add_dense(SparseT const &sval)
    {
        auto ii(_s2d.find(sval));
        if (ii == _s2d.end()) {
            auto densei(dense_extent());
            _s2d.insert(std::pair<SparseT,DenseT>(sval, densei));
            _d2s.push_back(sval);
            return densei;
        } else {
            return ii->second;
        }
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
enum class SparseTransform {
    ID,           // No transformation
    TO_DENSE,     // Convert sparse to dense indices
    TO_SPARSE,    // Convert dense to sparse indices
    ADD_DENSE     // Convert sparse to dense, adding to the SparseSet if it's not already there.
};

// =======================================================

template<class AccumT,class SparseSetT>
class SparseTransformAccum : public accum::Filter<AccumT>
{
    typedef accum::Filter<AccumT> super;
public:
    struct Data {
        SparseSetT * const sparse_set;
        SparseTransform transform;

        Data(SparseSetT *_sparse_set, SparseTransform _transform) :
            sparse_set(_sparse_set), transform(_transform) {}
    };

private:
    std::vector<Data> data;

public:
    /** @param _transform
        SparseTransform::TO_DENSE
        SparseTransform::TO_SPARSE
    @param sparse_sets
        0 if that dimension is not to be transformed
    */
    SparseTransformAccum(
        SparseTransform _transform,
        std::array<SparseSetT *, (size_t)super::rank> const sparse_sets,
        AccumT &&_sub);

    /** Merge shape we're given with shape from transforms. */
    void set_shape(std::array<long, super::rank> shape)
    {
        for (int i=0; i<super::rank; ++i) {
            // Use the shape we were given, if no transform for this dimension
            if (!data[i].sparse_set) continue;

            switch(data[i].transform) {
                case SparseTransform::ADD_DENSE:
                    // set_shape() is not very useful with ADD_DENSE;
                    // you will have to call it again after all
                    // items have been added.
                    // (fall through)
                case SparseTransform::TO_DENSE:
                    shape[i] = data[i].sparse_set->dense_extent();
                    break;
                case SparseTransform::TO_SPARSE:
                    shape[i] = data[i].sparse_set->sparse_extent();
                    break;
            }
        }
        super::sub.set_shape(shape);
    }

    void add(std::array<typename super::index_type,super::rank> index, typename super::val_type const &val)
    {
        for (int i=0; i<super::rank; ++i) {
            // Use the index we were given, if no transform for this dimension
            if (!data[i].sparse_set) continue;

            switch(data[i].transform) {
                case SparseTransform::ADD_DENSE:
                    index[i] = data[i].sparse_set->add_dense(index[i]);
                    break;
                case SparseTransform::TO_DENSE:
                    index[i] = data[i].sparse_set->to_dense(index[i]);
                    break;
                case SparseTransform::TO_SPARSE:
                    index[i] = data[i].sparse_set->to_sparse(index[i]);
                    break;
            }
        }
        super::sub.add(index, val);
    }
};

// ----------------------------------------------------------------
template<class AccumT,class SparseSetT>
SparseTransformAccum<AccumT,SparseSetT>::SparseTransformAccum(
    SparseTransform transform,
    std::array<SparseSetT *, (size_t)super::rank> const sparse_sets,    // 0 if we don't want to densify/sparsify that dimension
    AccumT &&_sub)
    : super(std::move(_sub))
{
    for (int i=0; i<super::rank; ++i) {
        data.push_back(Data(sparse_sets[i],
            sparse_sets[i] ? transform : SparseTransform::ID));
    }
}

namespace accum {

#define SparseTransformAccumT SparseTransformAccum<AccumT,SparseSetT>
template<class AccumT,class SparseSetT>
inline SparseTransformAccumT sparse_transform(
        SparseTransform transform,
        std::array<SparseSetT *, AccumT::rank> const &sparse_sets,
        AccumT &&sub)
{
    return SparseTransformAccumT(transform, sparse_sets, std::move(sub));
}
#undef SparseTransformAccumT

}    // namespace accum
// ----------------------------------------------------------------
template<class AccumT, class SrcT, class SparseT, class DenseT>
void sparse_copy(AccumT &ret, SrcT const &src,
    SparseTransform direction,
    std::array<spsparse::SparseSet<SparseT, DenseT> *, AccumT::rank> dims,
    bool set_shape=true);

template<class AccumT, class SrcT, class SparseT, class DenseT>
void sparse_copy(AccumT &ret, SrcT const &src,
    SparseTransform direction,
    std::array<spsparse::SparseSet<SparseT, DenseT> *, AccumT::rank> dims,
    bool set_shape)
{
    auto accum(sparse_transform_accum(&ret, direction, dims));
    spcopy(accum, src, set_shape);
}


}   // Namespace
