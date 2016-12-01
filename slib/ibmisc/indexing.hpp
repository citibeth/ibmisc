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


/** Fields describing one dimension */
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

/** Base class for Indexing<..> template: most of the core functionality,
without to_tuple() and to_index() methods.  Functions that can rely on
Indexing don't need to be templates themselves. */
class Indexing
{
protected:
    std::vector<IndexingData> data;

    /** Index IDs sorted by descending stride. {0,1,...} for row-major,
    reversed for col-major.  Can serve as a permutation (see
    permute.hpp) if you want the indices listed by descending
    stride.
    (NOTE: This is the reverse of Blitz stor.ordering()) */
    std::vector<int> indices;

    /** Computes the stride for each dimension. */
    void make_strides();

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
        std::vector<IndexingData> &&_data,
        std::vector<int> &&_indices);

    /** Number of dimensions in this indexing */
    size_t rank() const { return data.size(); }

    /** Total number of elements implied by this indexing.
    (Multiply extent of each dimension) */
    long extent() const;

    /** @return Information on the ix'th dimension */
    IndexingData const &operator[](size_t ix) const
        { return data[ix]; }

    /** Read/Write this indexing spec to a NetCDF file. */
    void ncio(NcIO &ncio, std::string const &vname);

    /** Allocates a blitz::Array according to the indexing described here. */
    template<class ValT, int RANK>
    blitz::Array<ValT,RANK> make_blitz();

    /** Creates a blitz::Array on existing memory, according to our indexing. */
    template<class ValT, int RANK>
    blitz::Array<ValT,RANK> to_blitz(ValT *data);

    template<class TupleT>
    long tuple_to_index(TupleT const *tuple) const
    {
        IndexT ix = 0;
        for (int k=0; k<rank(); ++k)
            ix += (tuple[k]-base[k]) * strides[k];
        return ix;
    }

    template<class TupleT>
    long tuple_to_index(std::vector<TupleT> const &tuple) const
        { return tuple_to_index(&tuple[0]); }

    template<class TupleT>
    long tuple_to_index(std::array<TupleT> const &tuple) const
        { return tuple_to_index(&tuple[0]); }

    template<class TupleT>
    void index_to_tuple(TupleT *tuple, long ix) const
    {
        for (int d=0; d< rank()-1; ++d) {       // indices by descending stride
            int const k = indices[d];
            TupleT tuple_k = ix / strides[k];
            ix -= tuple_k*strides[k];
            tuple[k] = tuple_k + base[k];
        }
        tuple[indices[rank()-1]] = ix;
    }


    template<class TupleT, int RANK>
    std::array<TupleT,RANK> index_to_tuple(IndexT ix) const
    {
        std::array<TupleT, RANK> ret;
        index_to_tuple(&ret[0], ix);
        return ret;
    }
};


/** By convention, indices are always in "alphabetical" order,
    regardless of the underlying storage.  This may be
    Fortran-order or C-order, depending.  Thus:

    ModelE: im,jm,ihc
    Ice model: ix,iy
*/
template<class TupleT>
Indexing make_indexing(
    std::vector<std::string> const &_name,
    std::vector<TupleT> const &_base,
    std::vector<TupleT> const &_extent,
    std::vector<int> &&_indices)
{
    std::vector<IndexingData> data;
    for (size_t i=0; i<_base.size(); ++i)
        data.push_back(IndexingData(_name[i], _base[i], _extent[i]));
    Indexing ret(std::move(data), std::move(indices));
    make_strides();
    return ret;
}




/** Adds the dimensions specified by indexing to a NcDimSpec.
@param permutation Permutation to apply to dims from indexing.
       If none given, then indexing.indices will be used
       (resulting in dimensions in decreasing stride order) */
NcDimSpec &append(NcDimSpec &dim_spec, Indexing &indexing,
    std::vector<int> const &_permutation = {})
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
blitz::Array<ValT,RANK> Indexing::make_blitz()
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
blitz::Array<ValT,RANK> Indexing::to_blitz(ValT *data)
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
// ==============================================================
struct DomainData {
    long const low;    // First "included" element in each index
    long const high;   // First "excluded" element in each index

    DomainData(long _low, long _high) : low(_low), high(_high) {}
};

/** Defines the boundaries of an MPI domain.  A simple hypercube in n-D space... */
class Domain {
public:
    std::vector<DomainData> data;

    Domain() {}

    template<class TupleT>
    Domain(std::vector<TupleT> const &_low, std::vector<TupleT> const &_high);

    int rank() const { return data.size(); }

    void ncio(NcIO &ncio, std::string const &vname);

    template<class TupleT>
    bool in_domain(TupleT const *tuple) const;

    template<class TupleT, int RANK>
    bool in_domain(std::array<TupleT, RANK> const &tuple) const
        { return in_domain(&tuple[0]); }

};

template<class TupleT>
Domain::Domain(std::vector<TupleT> const &_low, std::vector<TupleT> const &_high)
{
    for (size_t i=0; i<_low.size(); ++i)
        data.push_back(DomainData(_low[i], _high[i]));
}

template<class TupleT>
bool Domain::in_domain(TupleT const *tuple) const
{
    for (int k=0; k<rank(); ++k) {
        auto &datak(data[k]);
        if ((tuple[k] < datak.low) || (tuple[k] >= datak.high))
            return false;
    }
    return true;
}

extern bool in_domain(
    Domain const &domain,
    BaseIndexing const &indexing,
    IndexT ix);

// ============================================



} // Namespace
#endif
