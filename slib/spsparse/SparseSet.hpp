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
#include <ibmisc/netcdf.hpp>
#include <ibmisc/permutation.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/blitz.hpp>

namespace spsparse {

/** Translates between a sparse set (say, the set of indices used in a
SparseVector) and a dense set numbered [0...n) */
template<class SparseT, class DenseT>
class SparseSet {
public:
    SparseT _sparse_extent;
    std::unordered_map<SparseT, DenseT> _s2d;
    std::vector<SparseT> _d2s;

public:
    typedef SparseT sparse_type;
    typedef DenseT dense_type;


    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    SparseSet() : _sparse_extent(-1) {}

    void clear();

    bool in_sparse(SparseT const &sparse_ix) const
        { return _s2d.find(sparse_ix) != _s2d.end(); }

    bool in_dense(DenseT dense_ix) const
        { return (dense_ix >= 0 && dense_ix < dense_extent()); }

    SparseT sparse_extent() const
        { return _sparse_extent; }

    void set_sparse_extent(SparseT extent);

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
    {
        auto ii(_s2d.find(sval));
        if (ii == _s2d.end()) {
            std::stringstream sval_s;
            sval_s << sval;
            (*ibmisc::ibmisc_error)(-1,
                "Sparse value %s not found in SparseSet", sval_s.str().c_str());
        }
        return ii->second;
    }

    SparseT to_sparse(DenseT const &dval) const
    {
        if (dval < 0 || dval >= _d2s.size()) {
            (*ibmisc::ibmisc_error)(-1,
                "Value %ld is out of range (0, %ld)", (long)dval, _d2s.size());
        }
        return _d2s[dval];
    }


};

template<class SparseT, class DenseT>
void SparseSet<SparseT, DenseT>::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    if (ncio.rw == 'r') (*ibmisc::ibmisc_error)(-1,
        "SparseSet::ncio() currently does not support reading.");

    get_or_add_dim(ncio, vname + ".sparse_extent", _sparse_extent);
    // Set up dimensions
    auto dims(ibmisc::get_or_add_dims(ncio,
        {vname + ".dense_extent"},
        {dense_extent()}));

    ibmisc::ncio_vector(ncio, _d2s, false, vname, ibmisc::get_nc_type<SparseT>(), dims);
}

template<class SparseT, class DenseT>
void SparseSet<SparseT, DenseT>::clear() {
    _sparse_extent = -1;
    _s2d.clear();
    _d2s.clear();
}

template<class SparseT, class DenseT>
void SparseSet<SparseT, DenseT>::set_sparse_extent(SparseT extent)
{
    std::stringstream buf;
    if (_sparse_extent != -1 && _sparse_extent != extent) {
        buf << "Cannot change sparse_extent from " << _sparse_extent << " to " << extent;
        (*ibmisc::ibmisc_error)(-1, "%s", buf.str().c_str());
    }
    _sparse_extent = extent;
}
// ==================================================================
/** Creates a SparseSet instances that is intentionally an
    identity function. */
template<class SparseSetT>
SparseSetT id_sparse_set(typename SparseSetT::dense_type n)
{
    SparseSetT sparse_set;
    for (typename SparseSetT::dense_type i=0; i<n; ++i) sparse_set.add_dense(i);
    sparse_set.set_sparse_extent(sparse_set.dense_extent());
    return sparse_set;
}

// -----------------------------------------------------------
// This is in priority order to process; see process_order below.
enum class SparsifyTransform {
    TO_DENSE_IGNORE_MISSING,
    TO_DENSE,     // Convert sparse to dense indices
    TO_SPARSE,    // Convert dense to sparse indices
    ADD_DENSE,    // Convert sparse to dense, adding to the SparseSet if it's not already there.
    ID           // No transformation

};


// =======================================================

namespace accum {

// -----------------------------------------------------------
/** Accumulator to create a SparseSet of the indices used */
template<class SparseSetT, class ValueT, int RANK>
class SparseSetAccum {
    std::array<SparseSetT *, RANK> dims;
public:
    SparseSetAccum(std::array<SparseSetT *, RANK> _dims)
        : dims(_dims) {}

    void add(std::array<typename SparseSetT::sparse_type,RANK> const &index, ValueT const &val)
    {
        for (int i=0; i<RANK; ++i) if (dims[i]) dims[i]->add_dense(index[i]);
    }

};

template<class SparseSetT, class ValueT, int RANK>
SparseSetAccum<SparseSetT, ValueT, RANK> sparse_set(std::array<SparseSetT *, RANK> _dims)
    { return SparseSetAccum<SparseSetT,ValueT,RANK>(_dims); }

// -----------------------------------------------------------
#define SPARSIFY_TPARAMS class AccumT, class SparseSetT, class InIndexT
#define SparsifyT Sparsify<AccumT,SparseSetT,InIndexT>

template<SPARSIFY_TPARAMS>
class Sparsify : public Filter<AccumT>
{
    typedef Filter<AccumT> super;
    typedef typename super::index_type out_index_type;
public:
    typedef InIndexT index_type;    // Override super
public:
    struct Data {
        SparseSetT * const sparse_set;
        SparsifyTransform transform;

        Data(SparseSetT *_sparse_set, SparsifyTransform _transform) :
            sparse_set(_sparse_set), transform(_transform) {}

        long extent() const {
            switch(transform) {
                case SparsifyTransform::ID :
                    return -1;    // sparse_set == NULL
                case SparsifyTransform::ADD_DENSE :
                case SparsifyTransform::TO_DENSE_IGNORE_MISSING :
                case SparsifyTransform::TO_DENSE :
                    return sparse_set->dense_extent();
                case SparsifyTransform::TO_SPARSE :
                    return sparse_set->sparse_extent();
            }
        }
    };
private:

    std::vector<Data> data;

    // The order in which we will process dimensions.
    // This way, we can process dimensions that don't change
    // state (eg TO_DENSE_IGNORE_MISSING) before dimensions
    // that do change state.
    std::array<int, (size_t)super::rank> process_order;

public:
    Data &dim(int ix)
        { return data[ix]; }

    /** @param _transform.  Either a single item (used for all dimensions), or one per dimension.
        SparsifyTransform::TO_DENSE
        SparsifyTransform::TO_SPARSE
    @param sparse_sets
        0 if that dimension is not to be transformed
    */
    Sparsify(
        std::vector<SparsifyTransform> const &_transform,
        std::array<SparseSetT *, (size_t)super::rank> const sparse_sets,
        AccumT &&_sub);

    /** Merge shape we're given with shape from transforms. */
    void set_shape(std::array<long, super::rank> shape)
    {
        for (int i=0; i<super::rank; ++i) {
            // Use the shape we were given, if no transform for this dimension
            if (!data[i].sparse_set) continue;

            switch(data[i].transform) {
                case SparsifyTransform::ADD_DENSE:
                    // set_shape() is not very useful with ADD_DENSE;
                    // you will have to call it again after all
                    // items have been added.
                    // (fall through)
                case SparsifyTransform::TO_DENSE_IGNORE_MISSING:
                case SparsifyTransform::TO_DENSE:
                    shape[i] = data[i].sparse_set->dense_extent();
                    break;
                case SparsifyTransform::TO_SPARSE:
                    shape[i] = data[i].sparse_set->sparse_extent();
                    break;
            }
        }
        super::sub.set_shape(shape);
    }

    void add(std::array<index_type,super::rank> index, typename super::val_type const &val)
    {
        std::array<out_index_type, super::rank> index2;

        // Iterate in reverse order to accommodate the common case of
        // {ADD_DENSE, TO_DENSE_IGNORE_MISSING}.  This allows one to create
        // a matrix based on a subset of the input, adding only the output
        // grid cells that are needed.
        for (int ix=0; ix<super::rank; ++ix) {
            int const i = process_order[ix];

            switch(data[i].transform) {
                case SparsifyTransform::ID:
                    // Use the index given, no transform
                    index2[i] = index[i];
                    break;
                case SparsifyTransform::ADD_DENSE:
                    index2[i] = data[i].sparse_set->add_dense(index[i]);
                    break;
                case SparsifyTransform::TO_DENSE_IGNORE_MISSING: {
                    auto &sset(*data[i].sparse_set);
                    auto ii(sset._s2d.find(index[i]));
                    if (ii == sset._s2d.end()) return;    // An index was missing; ignore this element in the sparse matrix
                    index2[i] = ii->second;
                } break;
                case SparsifyTransform::TO_DENSE:
                    index2[i] = data[i].sparse_set->to_dense(index[i]);
                    break;
                case SparsifyTransform::TO_SPARSE:
                    index2[i] = data[i].sparse_set->to_sparse(index[i]);
                    break;
            }
        }
        super::sub.add(index2, val);
    }
};

// ----------------------------------------------------------------
template<SPARSIFY_TPARAMS>
SparsifyT::Sparsify(
    std::vector<SparsifyTransform> const &transforms,
    std::array<SparseSetT *, (size_t)super::rank> const sparse_sets,    // 0 if we don't want to densify/sparsify that dimension
    AccumT &&_sub)
    : super(std::move(_sub))
{
    if ((transforms.size() < 1) || (transforms.size() > 1 && transforms.size() != super::rank))
        (*ibmisc::ibmisc_error)(-1, "Must give 1 transform type per dimension, or 1 for all dimensions");

    // Tell what to do with each rank
    for (int i=0; i<super::rank; ++i) {
        SparsifyTransform transform = (i >= transforms.size() ? transforms[0] : transforms[i]);
        data.push_back(Data(sparse_sets[i],
            sparse_sets[i] ? transform : SparsifyTransform::ID));
    }

    // Determine processing order of the ranks
    // Iterate in reverse order to accommodate the common case of
    // {ADD_DENSE, TO_DENSE_IGNORE_MISSING}.  This allows one to create
    // a matrix based on a subset of the input, adding only the output
    // grid cells that are needed.
    std::vector<int> _process_order;
    ibmisc::sorted_permutation(data.begin(), data.end(), _process_order,
        [&](Data const &a, Data const &b)
            { return a.transform < b.transform; }
    );
    std::copy_n(_process_order.begin(), super::rank, process_order.begin());
}

template<class X>
struct in_index_type {
    typedef X x;
};

template<SPARSIFY_TPARAMS>
inline SparsifyT sparsify(
        std::vector<SparsifyTransform> const &transforms,
        in_index_type<InIndexT> const dummy,
        std::array<SparseSetT *, AccumT::rank> const &sparse_sets,
        AccumT &&sub)
{
    return SparsifyT(transforms, sparse_sets, std::move(sub));
}
#undef SparsifyT
#undef SPARSIFY_TPARAMS

// -------------------------------------------------------------
#define SPARSIFY_TPARAMS class AccumT, class SparseSetT
#define ToDenseT Sparsify<AccumT,SparseSetT,typename SparseSetT::sparse_type>
#define ToSparseT Sparsify<AccumT,SparseSetT,typename SparseSetT::dense_type>

template<SPARSIFY_TPARAMS>
inline ToDenseT to_dense(
        std::array<SparseSetT *, AccumT::rank> const &dims,
        AccumT &&sub)
{ return ToDenseT({SparsifyTransform::TO_DENSE}, dims, std::move(sub)); }

template<SPARSIFY_TPARAMS>
inline ToDenseT to_dense_ignore_missing(
        std::array<SparseSetT *, AccumT::rank> const &dims,
        AccumT &&sub)
{ return ToDenseT({SparsifyTransform::TO_DENSE_IGNORE_MISSING}, dims, std::move(sub)); }

template<SPARSIFY_TPARAMS>
inline ToDenseT add_dense(
        std::array<SparseSetT *, AccumT::rank> const &dims,
        AccumT &&sub)
{ return ToDenseT({SparsifyTransform::ADD_DENSE}, dims, std::move(sub)); }

template<SPARSIFY_TPARAMS>
inline ToSparseT to_sparse(
        std::array<SparseSetT *, AccumT::rank> const &dims,
        AccumT &&sub)
{ return ToSparseT({SparsifyTransform::TO_SPARSE}, dims, std::move(sub)); }

#undef ToSparseT
#undef ToDenseT
#undef SPARSIFY_TPARAMS



// ----------------------------------------------------------------

}    // namespace accum


}   // Namespace
