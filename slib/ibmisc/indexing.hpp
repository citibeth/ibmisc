#ifndef IBMISC_INDEXING
#define IBMISC_INDEXING

#include <ibmisc/netcdf.hpp>

namespace ibmisc {


template<class TupleT, class IndexT>
class Indexing
{
public:
	std::vector<TupleT> base;	// First element in each index
	std::vector<TupleT> extent;	// Extent (# elements) of each index
	std::vector<int> indices;	// Index IDs sorted by descending stride. {0,1,...} for row-major, reversed for row-major

	// Derived fields...
	std::vector<IndexT> strides;


	size_t rank() const { return extent.size(); }

protected:
	std::vector<IndexT> make_strides()
	{
		std::vector<IndexT> ret(rank());
		ret[indices[rank()-1]] = 1;
		for (int d=rank()-2; d>=0; --d) {
			ret[indices[d]] = ret[indices[d+1]] * extent[indices[d+1]];
		}
		return ret;
	}

public:
	Indexing() {}

	Indexing(
		std::vector<TupleT> &&_base,
		std::vector<TupleT> &&_extent,
		std::vector<int> &&_indices)
	: base(std::move(_base)),
		extent(std::move(_extent)),
		indices(std::move(_indices)),
		strides(make_strides())
	{}

	IndexT size() const
	{
		IndexT ret = extent[0];
		for (int k=1; k<rank(); ++k) ret *= extent[k];
		return ret;
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
		for (int d=0; d< rank()-1; ++d) {		// indices by descending stride
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


	void ncio(NcIO &ncio, netCDF::NcType ncTupleT, std::string const &vname);
};

template<class TupleT, class IndexT>
void Indexing<TupleT, IndexT>::ncio(
	NcIO &ncio,
	netCDF::NcType ncTupleT,
	std::string const &vname)
{
	auto info_v = get_or_add_var(ncio, vname, netCDF::ncInt64, {});
	get_or_put_att(info_v, ncio.rw, "base", ncTupleT, base);
	get_or_put_att(info_v, ncio.rw, "extent", ncTupleT, extent);
	get_or_put_att(info_v, ncio.rw, "indices", ncTupleT, indices);
}

// ----------------------------------------------------------------
/** Defines the boundaries of an MPI domain.  A simple hypercube in n-D space... */
template<class TupleT>
class Domain {
public:
	std::vector<TupleT> low;	// First "included" element in each index
	std::vector<TupleT> high;	// First "excluded" element in each index

	Domain() {}
	Domain(std::vector<TupleT> &&_low, std::vector<TupleT> &&_high)
		: low(std::move(_low)), high(std::move(_high)) {}

	int rank() const { return low.size(); }

	bool in_domain(TupleT *tuple) const
	{
		for (int k=0; k<rank; ++k) {
			if ((tuple[k] < low[k]) || (tuple[k] >= high[k])) return false;
		}
		return true;
	}

	template<int RANK>
	bool in_domain(std::array<TupleT, RANK> const &tuple) const
		{ return in_domain(&tuple[0]); }

	void ncio(NcIO &ncio, netCDF::NcType ncTupleT, std::string const &vname);
};

template<class TupleT>
void Domain<TupleT>::ncio(
	NcIO &ncio,
	netCDF::NcType ncTupleT,
	std::string const &vname)
{
	auto info_v = get_or_add_var(ncio, vname, netCDF::ncInt64, {});
	get_or_put_att(info_v, ncio.rw, "low", ncTupleT, low);
	get_or_put_att(info_v, ncio.rw, "high", ncTupleT, high);
}

template<class TupleT, class IndexT>
bool in_domain(
	Domain<TupleT> const *domain,
	Indexing<TupleT, IndexT> const *indexing,
	IndexT ix)
{
	TupleT tuple[domain.rank()];
	indexing.index_to_tuple(tuple, ix);
	return domain.in_domain(tuple);
}

// ============================================



} // Namespace
#endif
