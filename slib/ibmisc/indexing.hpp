#ifndef IBMISC_INDEXING
#define IBMISC_INDEXING

namespace ibmisc {

template<class TupleT, int RANK, class IndexT>
class Indexing {
public:

	virtual ~Indexing() {};
	virtual IndexT tuple_to_index(std::array<TupleT, RANK>) const = 0;
	virtual std::array<TupleT, RANK> index_to_tuple(IndexT index) const = 0;
};
// ---------------------------------------
template<class TupleT, int RANK, class IndexT>
class Indexing_ColMajor : public Indexing<TupleT, RANK, IndexT>
{
protected:
	std::array<TupleT, RANK> extent;
	std::array<IndexT, RANK> strides;

public:

	virtual ~Indexing_ColMajor() {}
	Indexing_ColMajor(std::array<TupleT, RANK> _extent) : extent(_extent)
	{
		strides[0] = 1;
		for (TupleT k=1; k<RANK; ++k) strides[k] = strides[k-1] * extent[k-1];
	}
	IndexT size()
	{
		IndexT ret = extent[0];
		for (int k=1; k<RANK; ++k) ret *= extent[k];
		return ret;
	}


	IndexT tuple_to_index(std::array<TupleT, RANK> tuple) const
	{
		IndexT ix = tuple[0];		// We know strides[0] == 1
		for (TupleT k=1; k<RANK; ++k) ix += tuple[k]*strides[k];
		return ix;
	}

	std::array<TupleT, RANK> index_to_tuple(IndexT ix) const
	{
		std::array<TupleT, RANK> tuple;
		for (TupleT k=RANK-1; k>0; --k) {
			tuple[k] = ix / strides[k];
			ix -= tuple[k]*strides[k];
		}
		tuple[0] = ix;
		return tuple;
	}
};
// ---------------------------------------
template<class TupleT, int RANK, class IndexT>
class Indexing_RowMajor : public Indexing<TupleT, RANK, IndexT>
{
protected:
	std::array<TupleT, RANK> extent;
	std::array<IndexT, RANK> strides;

public:

	virtual ~Indexing_RowMajor() {}
	Indexing_RowMajor(std::array<TupleT, RANK> _extent) : extent(_extent)
	{
		strides[RANK-1] = 1;
		for (TupleT k=RANK-2; k>=0; --k) strides[k] = strides[k+1] * extent[k+1];
	}
	IndexT size()
	{
		IndexT ret = extent[0];
		for (int k=1; k<RANK; ++k) ret *= extent[k];
		return ret;
	}

	IndexT tuple_to_index(std::array<TupleT, RANK> tuple) const
	{
		IndexT ix = tuple[RANK-1];		// We know strides[RANK-1] == 1
		for (TupleT k=RANK-2; k>=0; --k) ix += tuple[k]*strides[k];
		return ix;
	}

	std::array<TupleT, RANK> index_to_tuple(IndexT ix) const
	{
		std::array<TupleT, RANK> tuple;
		for (TupleT k=0; k<RANK-1; ++k) {
			tuple[k] = ix / strides[k];
			ix -= tuple[k]*strides[k];
		}
		tuple[RANK-1] = ix;
		return tuple;
	}
};



} // Namespace
#endif
