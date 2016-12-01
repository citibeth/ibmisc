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

#ifndef IBMISC_INDEXING
#define IBMISC_INDEXING

#include <ibmisc/netcdf.hpp>

namespace ibmisc {


struct IndexingData {
    std::string const name;    // Name of this dimension
    long const base;           // First element in this dimension
    long const extent;         // Number of elements in this dimension
private:
    long _stride;         // (DERIVED) Stride (in #elements) for this dimension
public:
    long stride() { return _stride; }

    IndexingData(std::string const &_name, long _base, long _extent) :
        name(_name), base(_base), extent(_extent) {}
};

class IndexingBase
{
protected:
    std::vector<IndexingData> data;

    // Index IDs sorted by descending stride. {0,1,...} for row-major,
    // reversed for col-major.  Can serve as a permutation (see
    // permute.hpp) if you want the indices listed by descending
    // stride.
    std::vector<int> indices;

    /** Computes the stride for each dimension. */
    void make_strides();

public:
    IndexingBase() {}

    /** By convention, indices are always in "alphabetical" order,
        regardless of the underlying storage.  This may be
        Fortran-order or C-order, depending.  Thus:

        ModelE: im,jm,ihc
        Ice model: ix,iy
    */
    IndexingBase(
        std::vector<IndexingData> &&_data,
        std::vector<int> &&_indices);

    size_t rank() const { return data.size(); }

    long extent() const;

    IndexingData const &operator[](size_t ix) const
        { return data[ix]; }

    void ncio(NcIO &ncio, netCDF::NcType nclong, std::string const &vname);

    /** Allocates a blitz::Array according to the indexing described here. */
    template<class ValT, int RANK>
    blitz::Array<ValT,RANK> make_blitz();

    /** Creates a blitz::Array on existing memory, according to our indexing. */
    template<class ValT, int RANK>
    blitz::Array<ValT,RANK> to_blitz(ValT *data);
};

/** Adds the dimensions specified by indexing to a NcDimSpec. */
NcDimSpec &append(NcDimSpec &dim_spec, IndexingBase &indexing,
    std::vector<int> const &_permutation)
{
    // Default permutation puts largest stride first for NetCDF
    std::vector<int> const *permutation =
        (_permutation.size() != 0 ? &permutation : &indexing.indices);

    for (size_t i=0; i<indexing.rank(); ++i) {
        int dimi = (*permutation)[i];
        dim_spec.append(indexing[dimi].name, indexing[dimi].extent);
    }

    return dim_spec;
}


/** Allocates a blitz::Array according to the indexing described here. */
template<class ValT, int RANK>
blitz::Array<ValT,RANK> IndexingBase::make_blitz()
{
    if (rank() != RANK) (*ibmisc_error)(-1,
        "Rank mismatch: %d vs %d", rank(), RANK);

    blitz::TinyVector<int, RANK> _lbounds, _extent;
    blitz::GeneralArrayStorage<RANK> _stor;
    for (int i=0; i<RANK; ++i) {
        _lbounds[i] = base[i];
        _extent[i] = extent[i];
        _stor.ordering()[RANK-i] = indices[i];
    }
    return blitz::Array<ArrT,RANK>(_lbounds, _extent, _stor);
}

/** Creates a blitz::Array on existing memory, according to our indexing. */
template<class ValT, int RANK>
blitz::Array<ValT,RANK> IndexingBase::to_blitz(ValT *data)
{
    if (rank() != RANK) (*ibmisc_error)(-1,
        "Rank mismatch: %d vs %d", rank(), RANK);

    blitz::TinyVector<int, RANK> _shape, _stride;
    blitz::GeneralArrayStorage<RANK> _stor;
    for (int i=0; i<RANK; ++i) {
        _shape[i] = extent[i];
        _stride[i] = strides[i];
        _stor.base()[i] = base[i];
        // Ordering is not needed because we're using stride
        // stor.ordering()[i] = i;      // Fortran ordering, blitz++ manual p31
    }

    return blitz::Array<ArrT,RANK>(data, _shape, _stride,
        blitz::neverDeleteData, stor);
}


// ---------------------------------------------------------------------
template<class TupleT, class IndexT>
class Indexing : public IndexingBase
{
public:

public:
    Indexing() : IndexingBase()
    {}

    /** By convention, indices are always in "alphabetical" order,
        regardless of the underlying storage.  This may be
        Fortran-order or C-order, depending.  Thus:

        ModelE: im,jm,ihc
        Ice model: ix,iy
    */
    Indexing(
        std::vector<std::string> const &_name,
        std::vector<TupleT> const &_base,
        std::vector<TupleT> const &_extent,
        std::vector<int> &&_indices)
    : indices(std::move(_indices))
    {
        for (size_t i=0; i<_base.size(); ++i)
            data.push_back(IndexingData(_name[i], _base[i], _extent[i]));
        make_strides();
    }

    IndexT tuple_to_index(TupleT const *tuple) const
    {
        IndexT ix = 0;
        for (int k=0; k<rank(); ++k)
            ix += (tuple[k]-base[k]) * strides[k];
        return ix;
    }

    IndexT tuple_to_index(std::vector<TupleT> const &tuple) const
        { return tuple_to_index(&tuple[0]); }

    void index_to_tuple(TupleT *tuple, IndexT ix) const
    {
        for (int d=0; d< rank()-1; ++d) {       // indices by descending stride
            int const k = indices[d];
            TupleT tuple_k = ix / strides[k];
            ix -= tuple_k*strides[k];
            tuple[k] = tuple_k + base[k];
        }
        tuple[indices[rank()-1]] = ix;
    }

    template<int RANK>
    IndexT tuple_to_index(std::array<TupleT, RANK> const &tuple) const
        { return tuple_to_index(&tuple[0]); }

    template<int RANK>
    std::array<TupleT,RANK> index_to_tuple(IndexT ix) const
    {
        std::array<TupleT, RANK> ret;
        index_to_tuple(&ret[0], ix);
        return ret;
    }

};


// ----------------------------------------------------------------
/** Defines the boundaries of an MPI domain.  A simple hypercube in n-D space... */
class DomainBase {
public:
    std::vector<long> low;    // First "included" element in each index
    std::vector<long> high;   // First "excluded" element in each index

    Domain() {}

    int rank() const { return low.size(); }

    void ncio(NcIO &ncio, std::string const &vname);

};


template<class TupleT>
class Domain : DomainBase {
public:
    Domain() {}
    Domain(std::vector<TupleT> const &_low, std::vector<TupleT> const &_high)
    {
        for (size_t i=0; i<_low.size(); ++i) {
            low.push_back(_low[i]);
            high.push_back(_high[i]);
        }
    }

    bool in_domain(TupleT const *tuple) const
    {
        for (int k=0; k<rank(); ++k) {
            if ((tuple[k] < low[k]) || (tuple[k] >= high[k])) return false;
        }
        return true;
    }

    template<int RANK>
    bool in_domain(std::array<TupleT, RANK> const &tuple) const
        { return in_domain(&tuple[0]); }

};

extern bool in_domain(
    Domain const &domain,
    BaseIndexing const &indexing,
    IndexT ix);

// ============================================



} // Namespace
#endif
