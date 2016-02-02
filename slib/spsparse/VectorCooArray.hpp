#ifndef SPSPARSE_VECTOR_COOARRAY_HPP
#define SPSPARSE_VECTOR_COOARRAY_HPP

#include <spsparse/array.hpp>

namespace spsparse {

template<class IndexT, class ValT, int RANK>
class VectorCooArray
{
public:
	static const int rank = RANK;
	typedef IndexT index_type;
	typedef ValT val_type;
	typedef std::array<index_type, rank> indices_type;

	std::array<size_t, RANK> shape;		// Extent of each dimension
	void set_shape(std::array<size_t, RANK> const &_shape) { shape = _shape; }

protected:
	typedef VectorCooArray<IndexT, ValT, RANK> ThisVectorCooArrayT;
	std::array<std::vector<IndexT>, RANK> index_vecs;
	std::vector<ValT> val_vec;

	// OPTIONAL
	// If sort_order[0] != -1, this is the index of the beginning of
	// each different value of dimension sort_order[0] (i.e. beginning
	// of each row, if we're sorted row-major).
	bool dim_beginnings_set;
	std::vector<size_t> _dim_beginnings;

public:
	bool edit_mode;		// Are we in edit mode?
	std::array<int,RANK> sort_order;	// Non-negative elements if this is sorted

	VectorCooArray();

	VectorCooArray(std::array<size_t, RANK> const &_shape);

	std::unique_ptr<ThisVectorCooArrayT> new_blank() const
		{ return std::unique_ptr<ThisVectorCooArrayT>(new ThisVectorCooArrayT(shape)); }
	ThisVectorCooArrayT make_blank() const
		{ return ThisVectorCooArrayT(shape); }

	IndexT &index(int dim, size_t ix)
		{ return index_vecs[dim][ix]; }
	IndexT const &index(int dim, size_t ix) const
		{ return index_vecs[dim][ix]; }

	ValT &val(size_t ix)
		{ return val_vec[ix]; }
	ValT const &val(size_t ix) const
		{ return val_vec[ix]; }

	std::array<IndexT, RANK> index(int ix) const {
		std::array<IndexT, RANK> index_ret;
		for (int k=0; k<RANK; ++k) index_ret[k] = index(k, ix);
		return index_ret;
	}
	std::vector<IndexT> index_vec(int ix) const {
		std::vector<IndexT> ret;
		for (int k=0; k<RANK; ++k) ret.push_back(index(k, ix));
		return ret;
	}
	void set_index(int ix, std::array<IndexT, RANK> const &idx)
		{ for (int k=0; k<RANK; ++k) index(k, ix) = idx[k]; }



	blitz::Array<IndexT, 1> indices(int dim) const
		{ return ibmisc::to_blitz(index_vecs[dim]); }
	blitz::Array<ValT, 1> vals() const
		{ return ibmisc::to_blitz(val_vec); }

	// Move semantics
	VectorCooArray(VectorCooArray &&other);
	void operator=(ThisVectorCooArrayT &&other);

	// Copy semantics
	VectorCooArray(VectorCooArray const &other);
	void operator=(ThisVectorCooArrayT const &other);


	// -------------------------------------------------
	size_t size() const
		{ return val_vec.size(); }
	void clear();
	void reserve(size_t size);

	// -------------------------------------------------
	// --------------------------------------------------
	/** Standard STL-type iterator for iterating through a VectorSparseMatrix. */

	typedef CooIterator<const std::array<IndexT, RANK>, const IndexT, RANK, const ValT, const ThisVectorCooArrayT> const_iterator;
	typedef CooIterator<std::array<IndexT, RANK>, IndexT, RANK, ValT, ThisVectorCooArrayT> iterator;

	iterator begin(int ix = 0)
		{ return iterator(this, ix); }
	iterator end(int ix = 0)
		{ return iterator(this, size() + ix); }
	const_iterator cbegin(int ix = 0) const
		{ return const_iterator(this, ix); }
	const_iterator cend(int ix = 0) const
		{ return const_iterator(this, size() + ix); }
	const_iterator begin(int ix = 0) const
		{ return const_iterator(this, ix); }
	const_iterator end(int ix = 0) const
		{ return const_iterator(this, size() - ix); }

	// typedef DimIndexIter<IndexT, ValT, iterator> dim_iterator;
	typedef DimIndexIter<const IndexT, const ValT, const_iterator> const_dim_iterator;

	const_dim_iterator dim_iter(int dim, int ix) const
		{ return const_dim_iterator(dim, const_iterator(this, ix)); }
	const_dim_iterator dim_begin(int dim) const
		{ return dim_iter(dim, 0); }
	const_dim_iterator dim_end(int dim) const
		{ return dim_iter(dim, size()); }
	// -------------------------------------------------
	/** Goes in to add mode: legal to add more things to the vector. */
	void edit()
	{
		edit_mode = true;
		sort_order[0] = -1;
	}

	void add(std::array<IndexT, RANK> const index, ValT const val);
	void add_blitz(blitz::TinyVector<IndexT, RANK> const &index, ValT const val);

	/** Mark that this is now in sorted form. */
	void set_sorted(std::array<int,RANK> _sort_order)
	{
		sort_order = _sort_order;
		edit_mode = false;
	}

	// --------------------------------------------------
	// In-place algos
	void consolidate(
		std::array<int, RANK> const &_sort_order,
		DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
		bool handle_nan = false);

	void transpose(std::array<int, RANK> const &sort_order)
	{
		OverwriteAccum<iterator> overwrite(begin());
		spsparse::transpose(overwrite, *this, sort_order);
	}

	blitz::Array<ValT, RANK> to_dense();

	// Sets and returns this->_dim_beginnings
	std::vector<size_t> const &dim_beginnings() const;

	DimBeginningsXiter<ThisVectorCooArrayT> dim_beginnings_xiter() const;

	std::ostream &operator<<(std::ostream &out) const;
};



// --------------------------- Method Definitions
template<class IndexT, class ValT, int RANK>
VectorCooArray<IndexT, ValT, RANK>::
	VectorCooArray() : edit_mode(true), dim_beginnings_set(false), sort_order() {
		sort_order[0] = -1;
		for (int k=0; k<RANK; ++k) shape[k] = 0;	// User must set this later
	}

template<class IndexT, class ValT, int RANK>
VectorCooArray<IndexT, ValT, RANK>::
	VectorCooArray(std::array<size_t, RANK> const &_shape)
	: shape(_shape), edit_mode(true), dim_beginnings_set(false), sort_order() {
		sort_order[0] = -1;
	}

template<class IndexT, class ValT, int RANK>
VectorCooArray<IndexT, ValT, RANK>::
	VectorCooArray(VectorCooArray &&other) :
		shape(other.shape),
		index_vecs(std::move(other.index_vecs)),
		val_vec(std::move(other.val_vec)),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(std::move(other._dim_beginnings)),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

template<class IndexT, class ValT, int RANK>
	void VectorCooArray<IndexT, ValT, RANK>::operator=(ThisVectorCooArrayT &&other) {
		shape = other.shape;
		index_vecs = std::move(other.index_vecs);
		val_vec = std::move(other.val_vec);
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = std::move(other._dim_beginnings);
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

template<class IndexT, class ValT, int RANK>
VectorCooArray<IndexT, ValT, RANK>::
	VectorCooArray(VectorCooArray const &other) :
		shape(other.shape),
		index_vecs(other.index_vecs),
		val_vec(other.val_vec),
		dim_beginnings_set(other.dim_beginnings_set),
		_dim_beginnings(other._dim_beginnings),
		edit_mode(other.edit_mode),
		sort_order(other.sort_order) {}

template<class IndexT, class ValT, int RANK>
	void VectorCooArray<IndexT, ValT, RANK>::operator=(ThisVectorCooArrayT const &other) {
		shape = other.shape;
		index_vecs = other.index_vecs;
		val_vec = other.val_vec;
		dim_beginnings_set = other.dim_beginnings_set;
		_dim_beginnings = other._dim_beginnings;
		edit_mode = other.edit_mode;
		sort_order = other.sort_order;
	}

template<class IndexT, class ValT, int RANK>
void VectorCooArray<IndexT, ValT, RANK>::clear() {
		for (int k=0; k<RANK; ++k) index_vecs[k].clear();
		val_vec.clear();
		dim_beginnings_set = false;
		_dim_beginnings.clear();
		edit_mode = true;
		sort_order[0] = -1;
	}

template<class IndexT, class ValT, int RANK>
void VectorCooArray<IndexT, ValT, RANK>::reserve(size_t size) {
		for (int k=0; k<RANK; ++k) index_vecs[k].reserve(size);
		val_vec.reserve(size);
	}

// ------------------------------------------------------------------------
template<class IndexT, class ValT, int RANK>
	void VectorCooArray<IndexT, ValT, RANK>::add(std::array<IndexT, RANK> const index, ValT const val)
	{
		if (!edit_mode) {
			(*spsparse_error)(-1, "Must be in edit mode to use VectorCooArray::add()");
		}

		// Check bounds
		for (int i=0; i<RANK; ++i) {
			if (index[i] < 0 || index[i] >= shape[i]) {
				std::ostringstream buf;
				buf << "Sparse index out of bounds: index=(";
				for (int j=0; j<RANK; ++j) {
					buf << index[j];
					buf << " ";
				}
				buf << ") vs. shape=(";
				for (int j=0; j<RANK; ++j) {
					buf << shape[j];
					buf << " ";
				}
				buf << ")";
				(*spsparse_error)(-1, buf.str().c_str());
			}
		}

		for (int i=0; i<RANK; ++i) index_vecs[i].push_back(index[i]);
		val_vec.push_back(val);
	}
// ------------------------------------------------------------------------
template<class IndexT, class ValT, int RANK>
	void VectorCooArray<IndexT, ValT, RANK>::add_blitz(blitz::TinyVector<IndexT, RANK> const &index, ValT const val)
	{
		if (!edit_mode) {
			(*spsparse_error)(-1, "Must be in edit mode to use VectorCooArray::add()");
		}

		// Check bounds
		for (int i=0; i<RANK; ++i) {
			if (index[i] < 0 || index[i] >= shape[i]) {
				std::ostringstream buf;
				buf << "Sparse index out of bounds: index=(";
				for (int j=0; j<RANK; ++j) {
					buf << index[j];
					buf << " ";
				}
				buf << ") vs. shape=(";
				for (int j=0; j<RANK; ++j) {
					buf << shape[j];
					buf << " ";
				}
				buf << ")";
				(*spsparse_error)(-1, buf.str().c_str());
			}
		}

		for (int i=0; i<RANK; ++i) index_vecs[i].push_back(index[i]);
		val_vec.push_back(val);
	}
// ---------------------------------------------------------------

template<class IndexT, class ValT, int RANK>
void VectorCooArray<IndexT, ValT, RANK>::consolidate(
		std::array<int, RANK> const &_sort_order,
		DuplicatePolicy duplicate_policy,
		bool handle_nan)
	{
		// Do nothing if we're already properly consolidated
		if (this->sort_order == _sort_order && !edit_mode) return;

		ThisVectorCooArrayT ret(shape);
		spsparse::consolidate(ret, *this, _sort_order, duplicate_policy, handle_nan);
		*this = std::move(ret);
	}

template<class IndexT, class ValT, int RANK>
	blitz::Array<ValT, RANK> VectorCooArray<IndexT, ValT, RANK>::to_dense()
	{
		blitz::Array<ValT, RANK> ret(ibmisc::to_tiny<int,size_t,rank>(shape));
		ret = 0;
		DenseAccum<ThisVectorCooArrayT> accum(ret);
		copy(accum, *this);
		return ret;
	}

template<class IndexT, class ValT, int RANK>
	// Sets and returns this->_dim_beginnings
	std::vector<size_t> const &VectorCooArray<IndexT, ValT, RANK>::dim_beginnings() const
	{
		// See if we need to compute it; lazy eval
		if (!dim_beginnings_set) {
			// Const cast OK here for lazy eval implementation
			ThisVectorCooArrayT *vthis = const_cast<ThisVectorCooArrayT *>(this);
			vthis->_dim_beginnings = spsparse::dim_beginnings(*this);
			vthis->dim_beginnings_set = true;
		}
		return _dim_beginnings;
	}

template<class IndexT, class ValT, int RANK>
	DimBeginningsXiter <VectorCooArray<IndexT, ValT, RANK>> VectorCooArray<IndexT, ValT, RANK>::dim_beginnings_xiter() const
	{
		auto &db(dim_beginnings());
		int const index_dim = sort_order[0];
		int const val_dim = sort_order[1];
		return DimBeginningsXiter<ThisVectorCooArrayT>(this, index_dim, val_dim, db.begin(), db.end());
	}


// ---------------------------------------------------------------------------
template<class IndexT, class ValT, int RANK>
std::ostream &operator<<(std::ostream &os, spsparse::VectorCooArray<IndexT, ValT, RANK> const &A)
	{ return spsparse::_ostream_out_array(os, A); }

template<class IndexT, class ValT>
using VectorCooMatrix = VectorCooArray<IndexT, ValT, 2>;

template<class IndexT, class ValT>
using VectorCooVector = VectorCooArray<IndexT, ValT, 1>;


}	// Namespace
#endif	// Guard
