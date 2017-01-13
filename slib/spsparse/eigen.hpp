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

#include <ibmisc/iter.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/blitz.hpp>
#include <Eigen/SparseCore>

namespace spsparse {

/** @defgroup eigen eigen.hpp

@brief Linkages between the SpSparse and Eigen libraries

@{
*/

// --------------------------------------------------------------
/** Copies an Eigen SparseMatrix into a rank-2 Spsparse accumulator. */
template<class AccumulatorT, class ValueT>
extern void spcopy(
    AccumulatorT &ret,
    Eigen::SparseMatrix<ValueT> const &M,
    bool set_shape=true);

template<class AccumulatorT, class ValueT>
void spcopy(
    AccumulatorT &ret,
    Eigen::SparseMatrix<ValueT> const &M,
    bool set_shape=true)
{
    if (set_shape) ret.set_shape({M.rows(), M.cols()});

    // See here for note on iterating through Eigen::SparseMatrix
    // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1104
    for (int k=0; k<M.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<ValueT>::InnerIterator ii(M,k); ii; ++ii) {
        ret.add({ii.row(), ii.col()}, ii.value());
    }}
}
// --------------------------------------------------------------
template<class AccumulatorT, int _Rows>
inline void spcopy(
    AccumulatorT &ret,
    Eigen::Matrix<typename AccumulatorT::val_type, _Rows, 1> const &M,    // eg. Eigen::VectorXd
    bool set_shape=true);


/** Copy from dense Eigen column vectors */
template<class AccumulatorT, int _Rows>
inline void spcopy(
    AccumulatorT &ret,
    Eigen::Matrix<typename AccumulatorT::val_type, _Rows, 1> const &M,    // eg. Eigen::VectorXd
    bool set_shape)
{
    for (size_t i=0; i<M.rows(); ++i) {
        ret.add({(typename AccumulatorT::index_type)i}, M(i,0));
    }
}

// --------------------------------------------------------------
/** Sum the rows or columns of an Eigen SparseMatrix.
@param dimi 0: sum rows, 1: sum columns */
template<class ValueT>
inline blitz::Array<ValueT,1> sum(
    Eigen::SparseMatrix<ValueT> const &M, int dimi)
{
    // Get our weight vector (in dense coordinate space)
    blitz::Array<double,1> ret;
    auto accum1(spsparse::blitz_accum_new(&ret));
    auto accum2(permute_accum(&accum1, in_rank<2>(), {dimi}));
    spcopy(accum2, M, true);

    return ret;
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
            i, i, invert ? 1./diag(i) : diag(i)));
    }

    Eigen::SparseMatrix<ValueT> M(diag.extent(0), diag.extent(0));
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}

// --------------------------------------------------------------


template<class ValueT>
inline Eigen::SparseMatrix<ValueT> weight_matrix(blitz::Array<ValueT,1> weights)
    { return diag_matrix(weights, false); }


template<class ValueT>
inline Eigen::SparseMatrix<ValueT> scale_matrix(blitz::Array<ValueT,1> weights)
    { return diag_matrix(weights, true); }
// --------------------------------------------------------------
template<class ValueT>
inline Eigen::SparseMatrix<ValueT> weight_matrix(Eigen::SparseMatrix<ValueT> const &M, int dimi)
    { return diag_matrix(sum(M, dimi), false); }

template<class ValueT>
inline Eigen::SparseMatrix<ValueT> scale_matrix(Eigen::SparseMatrix<ValueT> const &M, int dimi)
    { return diag_matrix(sum(M, dimi), true); }

// --------------------------------------------------------------
/** An N-dimensional generalization of Eigen::Triplet */
template<class IndexT, class ValueT, int RANK>
class Tuple {
    std::array<IndexT, RANK> _index;
    ValueT _value;
public:
    IndexT &index(int i)
        { return i; }
    std::array<IndexT, RANK> &index()
        { return _index; }

    ValueT &value()
        { return _value; }

    Tuple(std::array<IndexT, RANK> const &index,
        ValueT const &value)
        : _index(index), _value(value) {}

    // ----- For Eigen::SparseMatrix::setFromTriplets()
    ValueT row() const
        { return index(0); }
    ValueT col() const
        { return index(1); }
};


/** Serves as accumulator and iterable storage */
template<class IndexT, class ValueT, int RANK>
class TupleVector : public std::vector<Tuple<IndexT,ValueT,RANK>>
{
public:
    typedef std::vector<Tuple<IndexT,ValueT,RANK>> super;
    static const int rank = RANK;
    typedef IndexT index_type;
    typedef ValueT value_type;

    
    std::array<long, rank> shape;

public:
    TupleVector() {}

    TupleVector(std::array<long,RANK> _shape)
        { set_shape(_shape); }

    // So this can serve as a Spsparse Accumulator
    void set_shape(std::array<long, rank> _shape)
        { shape = _shape; }

    void add(std::array<index_type,rank> const &index, ValueT const &value)
    {
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
                (*ibmisc::ibmisc_error)(-1, buf.str().c_str());
            }
        }

        super::push_back(Tuple<IndexT,ValueT,RANK>(index, value));
    }

    // see: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1370
    Eigen::SparseMatrix<ValueT,0,IndexT> to_eigen()
    {
        Eigen::SparseMatrix<ValueT,0,IndexT> M(shape[0], shape[1]);
        M.setFromTriplets(super::begin(), super::end());
        return M;
    }
};

template<class IndexT, class ValueT>
using TripletVector = TupleVector<IndexT,ValueT,2>;

// ---------------------------------------------

#define ARGS _Scalar,_Options,_StorageIndex
template<class _Scalar, int _Options, class _StorageIndex>
class EigenSparseMatrixIterator :
    public ibmisc::forward_iterator<
        EigenSparseMatrixIterator<ARGS>,
        EigenSparseMatrixIterator<ARGS>>
{
    typedef EigenSparseMatrixIterator<ARGS> IteratorT;
public:
    typedef _StorageIndex index_type;
    typedef _Scalar value_type;
    static const int RANK = 2;

    Eigen::SparseMatrix<ARGS> const &M;
    int k;
    typename Eigen::SparseMatrix<ARGS>::InnerIterator ii;
    
    EigenSparseMatrixIterator(Eigen::SparseMatrix<ARGS> const &_M, int _k)
        : M(_M), k(_k),
        ii(typename Eigen::SparseMatrix<ARGS>::InnerIterator(M,k)) {}

    // ---------------------------------------------------------
    // http://www.cplusplus.com/reference/iterator/

    IteratorT &operator++() {    // Prefix ++
        ++ii;
        if (!ii) {
            ++k;
            ii = typename Eigen::SparseMatrix<ARGS>::InnerIterator(M,k);
        }
        return *this;
    }

    bool operator==(IteratorT &other)
        { return (k == other->k) && (ii == other->ii); }
//    bool operator!=(IteratorT &other)
//        { return !(this == other); }

    IteratorT &operator*()
        { return *this; }
//    IteratorT *operator->()
//        { return this; }

    // ---------- What you get once you dereference

    _Scalar const &value()
        { return ii.value(); }
    _StorageIndex const &row()
        { return ii.row(); }
    _StorageIndex const &col()
        { return ii.col(); }
    _StorageIndex const &index(int ix)
        { return ix == 0 ? row() : col(); }
    std::array<_StorageIndex,2> index()
        { return {row(), col()}; }
};

// -----------------------------------

template<class _Scalar, int _Options, class _StorageIndex>
EigenSparseMatrixIterator<ARGS> begin(Eigen::SparseMatrix<ARGS> const &M)
    { return EigenSparseMatrixIterator<ARGS>(M,0); }

template<class _Scalar, int _Options, class _StorageIndex>
EigenSparseMatrixIterator<ARGS> end(Eigen::SparseMatrix<ARGS> const &M)
    { return EigenSparseMatrixIterator<ARGS>(M,M.outerSize()); }

#undef ARGS
// -----------------------------------

/** @} */

}   // namespace

