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

#include <spsparse/SparseSet.hpp>
#include <spsparse/VectorCooArray.hpp>
#include <spsparse/accum.hpp>
#include <Eigen/SparseCore>

namespace spsparse {

/** @defgroup eigen eigen.hpp

@brief Linkages between the SpSparse and Eigen libraries: basically,
converting stuff to/from Eigen::SparseMatrix with dense index space.

@{
*/

// --------------------------------------------------------------
/** Converts an Eigen::SparseMatrix to a rank-2 SpSparse structure.
In the process, it translates for each dimension from the dense index
space used by Eigen to a sparse index space used by SpSparse.

@param ret OUT: Place to store the resulting SpSparse array.
@param M IN: The Eigen structure to convert.
@param dims IN: Dense-to-sparse index conversions for the two dimensions. */
template<class AccumulatorT, class EigenSparseMatrixT>
void eigenM_to_sparseM(
    AccumulatorT &ret,
    EigenSparseMatrixT const &M,
    std::array<SparseSet<
        typename AccumulatorT::index_type,
        typename EigenSparseMatrixT::Index> *, 2> const &dims);

template<class AccumulatorT, class EigenSparseMatrixT>
void eigenM_to_sparseM( // MM = Matrix->Matrix
    AccumulatorT &ret,
    EigenSparseMatrixT const &M,
    std::array<SparseSet<
        typename AccumulatorT::index_type,
        typename EigenSparseMatrixT::Index> *, 2> const &dims)
{
    ret.set_shape({dims[0]->sparse_extent(), dims[1]->sparse_extent()});
    for (int k=0; k<M.outerSize(); ++k) {
    for (typename EigenSparseMatrixT::InnerIterator ii(M,k); ii; ++ii) {
        ret.add({dims[0]->to_sparse(ii.row()), dims[1]->to_sparse(ii.col())}, ii.value());
    }}
}
// --------------------------------------------------------------

template<class SparseMatrixT>
struct _ES {
    typedef typename SparseMatrixT::index_type SparseIndexT;
    typedef typename SparseMatrixT::val_type ValT;
    typedef Eigen::SparseMatrix<ValT> EigenSparseMatrixT;
    typedef typename EigenSparseMatrixT::Index DenseIndexT;
    typedef SparseSet<SparseIndexT, DenseIndexT> SparseSetT;
    static const int rank = SparseMatrixT::rank;

    typedef spsparse::VectorCooArray<SparseIndexT, ValT, 1> SparseVectorT;
    typedef Eigen::Triplet<ValT> EigenTripletT;
};



/** A data sructure that acts like a normal SpSparse Accumulator.  BUT:
  1) When stuff is added to it, it also updates corresponding dimension maps (SparseSet).
  2) An Eigen::SparseMatrix may be extracted from it when construction is complete. */
template <class SparseMatrixT>
class SparseTriplets : public
    spsparse::MappedArray<SparseMatrixT, typename Eigen::SparseMatrix<typename SparseMatrixT::val_type>::Index>
{
public:
    typedef spsparse::MappedArray<SparseMatrixT, typename Eigen::SparseMatrix<typename SparseMatrixT::val_type>::Index> super;
    typedef typename super::SparseIndexT SparseIndexT;
    typedef typename super::ValT ValT;
    typedef typename super::SparseSetT SparseSetT;

    typedef spsparse::VectorCooArray<SparseIndexT, ValT, 1> SparseVectorT;
    typedef Eigen::SparseMatrix<ValT> EigenSparseMatrixT;
    typedef typename EigenSparseMatrixT::Index DenseIndexT;
    typedef Eigen::Triplet<ValT> EigenTripletT;

    /** @param dims Dimension maps for each dimension.  If one is
    preparing to multiply matrices, then each dimension map will be
    shared by at least two SparseTriplets object. */
    SparseTriplets(std::array<SparseSetT *,2> const &_dims) : super(_dims) {}

    /** Produces an Eigen::SparseMatrix from our internal SparseMatrixT.

    @param transpose Set to 'T' to make the result be a transpose of
        our internal M.  Use '.' for no transpose.
    @param invert Set to true to make the result be an element-wise
        multiplicative ineverse of M. */
    typename _ES<SparseMatrixT>::EigenSparseMatrixT to_eigen(char transpose='.', bool invert=false) const;


    /** Produces a diagonal Eigen::SparseMatrix S where S[i,i] is the sum of internal M over dimension 1-i.
    @param dimi Index of the dimension that should REMAIN in the scale matrix.
    @return Eigen::SparseMatrix of shape [len(dimi), len(dimi)] */
    typename _ES<SparseMatrixT>::EigenSparseMatrixT eigen_scale_matrix(int dimi) const;

};
// --------------------------------------------------------
template <class SparseMatrixT>
typename _ES<SparseMatrixT>::EigenSparseMatrixT SparseTriplets<SparseMatrixT>::to_eigen(char transpose, bool invert) const
{
    std::vector<EigenTripletT> triplets;

    for (auto ii=super::M.begin(); ii != super::M.end(); ++ii) {
        int ix0 = (transpose == 'T' ? 1 : 0);
        int ix1 = (transpose == 'T' ? 0 : 1);
        auto dense0 = super::dims[ix0]->to_dense(ii.index(ix0));
        auto dense1 = super::dims[ix1]->to_dense(ii.index(ix1));
        triplets.push_back(EigenTripletT(dense0, dense1, invert ? 1./ii.val() : ii.val()));
    }
    EigenSparseMatrixT ret(
        super::dims[transpose == 'T' ? 1 : 0]->dense_extent(),
        super::dims[transpose == 'T' ? 0 : 1]->dense_extent());
    ret.setFromTriplets(triplets.begin(), triplets.end());
    return ret;
}
// --------------------------------------------------------
template<class ValueT>
inline blitz::Array<ValueT> sum(
    Eigen::SparseMatrix<ValueT> const &M, int dimi)
{
    // Get our weight vector (in dense coordinate space)
    typename Eigen::SparseMatrix<ValueT>::Index dim_size = (dimi == 0 ? M.rows() : M.cols());
    blitz::Array<typename EigenSparseMatrixT::value_type,1> ret(dim_size);
    ret = 0;
    for (int k=0; k<M.outerSize(); ++k) {
    for (typename EigenSparseMatrixT::InnerIterator ii(M,k); ii; ++ii) {
        ret[dimi == 0 ? ii.row() : ii.col()] += ii.value();
    }}
}
// --------------------------------------------------------------
template<class ValueT>
Eigen::SparseMatrix<ValueT> diag_matrix(
    blitz::Array<ValueT,1> const &diag, bool invert);

template<class ValueT>
Eigen::SparseMatrix<ValueT> diag_matrix(
    blitz::Array<ValueT,1> const &diag, bool invert)
{
    int n = diag.extent(0);

    // Invert weight vector into Eigen-format scale matrix
    std::vector<Eigen::Triplet<ValueT>> triplets;
    for (int i=0; i < diag.extent(0); ++i) {
        triplets.push_back(Eigen::Triplet<ValueT>(
            i, i, invert ? 1./ii.val() : ii.val()));
    }

    EigenSparseMatrixT M(diag.extent(0), diag.extent(0));
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}

template<class ValueT>
Eigen::SparseMatrix<ValueT> diag_matrix(blitz::Array<ValueT,1> diag, bool invert)
inline weight_matrix(blitz::Array<ValueT,1> weights)
    { return diag_matrix(weights, false); }


template<class ValueT>
Eigen::SparseMatrix<ValueT> diag_matrix(blitz::Array<ValueT,1> diag, bool invert)
inline scale_matrix(blitz::Array<ValueT,1> weights)
    { return diag_matrix(weights, true); }
// --------------------------------------------------------------


/** @} */

}   // namespace

