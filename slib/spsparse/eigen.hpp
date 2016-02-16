#pragma once

#include <spsparse/SparseSet.hpp>
#include <spsparse/VectorCooArray.hpp>
#include <Eigen/SparseCore>

namespace spsparse {

template<class SparseMatrixT, class EigenSparseMatrixT>
void eigen_to_sparse(
	SparseMatrixT &ret,
	EigenSparseMatrixT const &M,
	std::array<SparseSet<
		typename SparseMatrixT::index_type,
		typename EigenSparseMatrixT::Index> *, 2> const &dims)
{
	ret.set_shape({dims[0]->sparse_extent(), dims[1]->sparse_extent()});
	for (int k=0; k<M.outerSize(); ++k) {
	for (typename EigenSparseMatrixT::InnerIterator ii(M,k); ii; ++ii) {
		ret.add({dims[0]->to_sparse(ii.row()), dims[1]->to_sparse(ii.col())}, ii.value());
	}}
}

template <class SparseMatrixT>
class SparseTriplets
{
public:
	typedef typename SparseMatrixT::index_type SparseIndexT;
	typedef typename SparseMatrixT::val_type ValT;
	typedef spsparse::VectorCooArray<SparseIndexT, ValT, 1> SparseVectorT;
	typedef Eigen::SparseMatrix<ValT> EigenSparseMatrixT;
	typedef typename EigenSparseMatrixT::Index DenseIndexT;
	typedef Eigen::Triplet<ValT> EigenTripletT;
	typedef SparseSet<SparseIndexT, DenseIndexT> SparseSetT;

	SparseMatrixT M;
	std::array<SparseSetT *, 2> const dims;

	SparseTriplets(std::array<SparseSetT *,2> const &_dims) : dims(_dims) {}

	void set_shape(std::array<size_t, 2> const &shape)
	{
		M.set_shape(shape);
		dims[0]->set_sparse_extent(shape[0]);
		dims[1]->set_sparse_extent(shape[1]);
	}

	void add(std::array<SparseIndexT, 2> const index, ValT const val)
	{
		M.add(index, val);
		dims[0]->add(index[0]);
		dims[1]->add(index[1]);
	}

	EigenSparseMatrixT to_eigen(char transpose='.', bool invert=false) const
	{
		std::vector<EigenTripletT> triplets;

		for (auto ii=M.begin(); ii != M.end(); ++ii) {
			int ix0 = (transpose == 'T' ? 1 : 0);
			int ix1 = (transpose == 'T' ? 0 : 1);
			auto dense0 = dims[ix0]->to_dense(ii.index(ix0));
			auto dense1 = dims[ix1]->to_dense(ii.index(ix1));
			triplets.push_back(EigenTripletT(dense0, dense1, invert ? 1./ii.val() : ii.val()));
		}
		EigenSparseMatrixT ret(
			dims[transpose == 'T' ? 1 : 0]->dense_extent(),
			dims[transpose == 'T' ? 0 : 1]->dense_extent());
		ret.setFromTriplets(triplets.begin(), triplets.end());
		return ret;
	}

	EigenSparseMatrixT eigen_scale_matrix(int dimi)
	{

		// Get our weight vector
		SparseVectorT weight({dims[dimi]->sparse_extent()});
		for (auto ii=M.begin(); ii != M.end(); ++ii) {
			weight.add({ii.index(dimi)}, ii.val());
		}
		weight.consolidate({0});

		// Convert to Eigen-format scale matrix (and convert indices to dense)
		std::vector<EigenTripletT> triplets;
		for (auto ii=weight.begin(); ii != weight.end(); ++ii) {
			auto dense = dims[dimi]->to_dense(ii.index(0));
			triplets.push_back(EigenTripletT(dense, dense, 1./ii.val()));
		}

		auto dense_extent(dims[dimi]->dense_extent());
		EigenSparseMatrixT scale(dense_extent, dense_extent);
		scale.setFromTriplets(triplets.begin(), triplets.end());
		return scale;

	}

};



template<class EigenSparseMatrixT>
EigenSparseMatrixT scale_matrix(EigenSparseMatrixT &M, int dimi)
{
	typedef VectorCooArray<typename EigenSparseMatrixT::Index, typename EigenSparseMatrixT::Scalar, 1> SparseVectorT;
	typedef Eigen::Triplet<typename EigenSparseMatrixT::Scalar, typename EigenSparseMatrixT::Index> EigenTripletT;

	// Get our weight vector (in dense coordinate space)
	typename EigenSparseMatrixT::Index dim_size = (dimi == 0 ? M.rows() : M.cols());
	SparseVectorT weight({dim_size});

	for (int k=0; k<M.outerSize(); ++k) {
	for (typename EigenSparseMatrixT::InnerIterator ii(M,k); ii; ++ii) {
		weight.add({dimi == 0 ? ii.row() : ii.col()}, ii.value());
	}}
	weight.consolidate({0});

	// Invert weight vector into Eigen-format scale matrix
	std::vector<EigenTripletT> triplets;
	for (auto ii=weight.begin(); ii != weight.end(); ++ii) {
		triplets.push_back(EigenTripletT(ii.index(0), ii.index(0), 1./ii.val()));
	}

	EigenSparseMatrixT scale(weight.shape[0], weight.shape[0]);
	scale.setFromTriplets(triplets.begin(), triplets.end());
	return scale;

}



}	// namespace

