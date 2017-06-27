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

#include <Eigen/SparseCore>
#include <ibmisc/iter.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/blitz.hpp>
#include <spsparse/SparseSet.hpp>

namespace spsparse {

/** @defgroup eigen eigen.hpp

@brief Linkages between the SpSparse and Eigen libraries

@{
*/

// --------------------------------------------------------------
#define ARGS _Scalar,_Options,_StorageIndex
/** Copies an Eigen SparseMatrix into a rank-2 Spsparse accumulator. */
template<class AccumT, class _Scalar, int _Options, class _StorageIndex>
extern void spcopy(
    AccumT &&ret,
    Eigen::SparseMatrix<ARGS> const &M,
    bool set_shape=true);

template<class AccumT, class _Scalar, int _Options, class _StorageIndex>
void spcopy(
    AccumT &&ret,
    Eigen::SparseMatrix<ARGS> const &M,
    bool set_shape=true)
{
    if (set_shape) ret.set_shape({M.rows(), M.cols()});

    // See here for note on iterating through Eigen::SparseMatrix
    // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1104
    for (int k=0; k<M.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<ARGS>::InnerIterator ii(M,k); ii; ++ii) {
        ret.add({ii.row(), ii.col()}, ii.value());
    }}
}
#undef ARGS
// --------------------------------------------------------------
template<class AccumT, int _Rows>
extern void spcopy(
    AccumT &&ret,
    Eigen::Matrix<typename AccumT::val_type, _Rows, 1> const &M,    // eg. Eigen::VectorXd
    bool set_shape=true);


/** Copy from dense Eigen column vectors */
template<class AccumT, int _Rows>
void spcopy(
    AccumT &&ret,
    Eigen::Matrix<typename AccumT::val_type, _Rows, 1> const &M,    // eg. Eigen::VectorXd
    bool set_shape)
{
    for (size_t i=0; i<M.rows(); ++i) {
        ret.add({(typename AccumT::index_type)i}, M(i,0));
    }
}

// --------------------------------------------------------------

template<class ValT>
Eigen::SparseMatrix<ValT> diag_matrix(
    blitz::Array<ValT,1> const &diag, char invert='+');

template<class ValT>
Eigen::SparseMatrix<ValT> diag_matrix(
    blitz::Array<ValT,1> const &diag, char invert='+')
{
    int n = diag.extent(0);

    // Invert weight vector into Eigen-format scale matrix
    std::vector<Eigen::Triplet<ValT>> triplets;
    for (int i=0; i < diag.extent(0); ++i) {
        if (diag(i) == 0) continue;
        triplets.push_back(Eigen::Triplet<ValT>(
            i, i, invert == '-' ? 1./diag(i) : diag(i)));
    }

    Eigen::SparseMatrix<ValT> M(diag.extent(0), diag.extent(0));
    M.setFromTriplets(triplets.begin(), triplets.end());

    return M;
}

// --------------------------------------------------------------
/** An N-dimensional generalization of Eigen::Triplet */
template<class IndexT, class ValT, int RANK>
class Tuple {
    std::array<IndexT, RANK> _index;
    ValT _value;
public:
    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & _index;
        ar & _value;
    }


    IndexT &index(int i)
        { return _index[i]; }
    IndexT const &index(int i) const
        { return _index[i]; }

    std::array<IndexT, RANK> &index()
        { return _index; }
    std::array<IndexT, RANK> const &index() const
        { return _index; }

    ValT &value()
        { return _value; }
    ValT const &value() const
        { return _value; }

    Tuple() {}    // WARNING: Uninitialized.  It's here for boost serialization

    Tuple(std::array<IndexT, RANK> const &index,
        ValT const &value)
        : _index(index), _value(value) {}

    // ----- For Eigen::SparseMatrix::setFromTriplets()
    ValT row() const
        { return index(0); }
    ValT col() const
        { return index(1); }

    bool operator<(Tuple<IndexT,ValT,RANK> const &other) const
        { return _index < other._index; }

};

/** Serves as accumulator and iterable storage */
template<class IndexT, class ValT, int RANK>
class TupleList : public std::vector<Tuple<IndexT,ValT,RANK>>
{
public:
    typedef std::vector<Tuple<IndexT,ValT,RANK>> super;
    static const int rank = RANK;
    typedef IndexT index_type;
    typedef ValT val_type;
    typedef TupleList base_array_type;

    
    std::array<long, rank> _shape;

public:
    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & *static_cast<super *>(this);
        ar & _shape;
    }

    base_array_type &base()
        { return *this; }
    base_array_type const &base() const
        { return *this; }

    struct iterator : public super::iterator {
        static const int rank = RANK;
        typedef IndexT index_type;
        typedef ValT val_type;

        iterator(typename super::iterator &&ii) : super::iterator(std::move(ii)) {}
    };
    iterator begin()
        { return iterator(super::begin()); }
    iterator end()
        { return iterator(super::end()); }
    // -------------------------------------------------
    struct const_iterator : public super::const_iterator {
        static const int rank = RANK;
        typedef IndexT index_type;
        typedef ValT val_type;

        const_iterator(typename super::const_iterator &&ii) : super::const_iterator(std::move(ii)) {}
    };
    const_iterator begin() const
        { return const_iterator(super::begin()); }
    const_iterator end() const
        { return const_iterator(super::end()); }

    // -------------------------------------------------


    TupleList()
        { _shape.fill(-1); }    // -1 means unlimited dimension

    TupleList(std::array<long,RANK> shape) : _shape(shape) {}

    // So this can serve as a Spsparse Accumulator
    void set_shape(std::array<long, rank> shape)
        { _shape = shape; }

    std::array<long,rank> const &shape() const
        { return _shape; }
    long shape(int i) const
        { return _shape[i]; }

    void add(std::array<index_type,rank> const &index, ValT const &value);
};

template<class IndexT, class ValT, int RANK>
void TupleList<IndexT,ValT,RANK>::add(std::array<index_type,rank> const &index, ValT const &value)
{
    // Check bounds
    for (int i=0; i<RANK; ++i) {
        if (_shape[i] >= 0 && (index[i] < 0 || index[i] >= _shape[i])) {
            std::ostringstream buf;
            buf << "Sparse index out of bounds: index=(";
            for (int j=0; j<RANK; ++j) {
                buf << index[j];
                buf << " ";
            }
            buf << ") vs. shape=(";
            for (int j=0; j<RANK; ++j) {
                buf << _shape[j];
                buf << " ";
            }
            buf << ")";
            (*ibmisc::ibmisc_error)(-1, buf.str().c_str());
        }
    }

    super::push_back(Tuple<IndexT,ValT,RANK>(index, value));
}

// --------------------------------------------------------

template<class ArrayT>
extern Eigen::SparseMatrix<typename ArrayT::val_type,0,typename ArrayT::index_type>
    to_eigen(ArrayT const &tuples);

template<class ArrayT>
Eigen::SparseMatrix<typename ArrayT::val_type,0,typename ArrayT::index_type>
    to_eigen(ArrayT const &tuples)
{
    Eigen::SparseMatrix<typename ArrayT::val_type,0,typename ArrayT::index_type> M(tuples.base().shape(0), tuples.base().shape(1));
    M.setFromTriplets(tuples.base().begin(), tuples.base().end());
    return M;
}


template<class AccumT, class IndexT, class ValT, int RANK>
extern void spcopy(AccumT &&ret, TupleList<IndexT,ValT,RANK> const &A, bool set_shape = true);

template<class AccumT, class IndexT, class ValT, int RANK>
void spcopy(AccumT &&ret, TupleList<IndexT,ValT,RANK> const &A, bool set_shape)
{
    if (set_shape) ret.set_shape(A.shape());
    for (auto ii=A.begin(); ii != A.end(); ++ii) {
        ret.add(ii->index(), ii->value());
    }
}
// ---------------------------------------------
// =========================================================

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
    typedef _Scalar val_type;
    static const int RANK = 2;

    Eigen::SparseMatrix<ARGS> const &M;
    int k;
    typename Eigen::SparseMatrix<ARGS>::InnerIterator ii;
    
    EigenSparseMatrixIterator(Eigen::SparseMatrix<ARGS> const &_M, int _k)
        : M(_M), k(_k),
        ii(typename Eigen::SparseMatrix<ARGS>::InnerIterator(M,k)) {}

    // ---------------------------------------------------------
    // http://www.cplusplus.com/reference/iterator/

    IteratorT &operator++();

    bool operator==(IteratorT const &other) const
    {
        if (k != other->k) return false;
        return (k == M.outerSize()) || (ii == other->ii);
    }
//    bool operator==(IteratorT const other) const
//        { return (k == other->k) && (ii == other->ii); }

    IteratorT &operator*()
        { return *this; }
    IteratorT const &operator*() const
        { return *this; }
//    IteratorT *operator->()
//        { return this; }

    // ---------- What you get once you dereference

    _Scalar const &value()
        { return ii.value(); }
    _StorageIndex row()
        { return ii.row(); }
    _StorageIndex col()
        { return ii.col(); }
    _StorageIndex index(int ix)
        { return ix == 0 ? row() : col(); }
    std::array<_StorageIndex,2> index()
        { return {row(), col()}; }
};

template<class _Scalar, int _Options, class _StorageIndex>
EigenSparseMatrixIterator<ARGS> &EigenSparseMatrixIterator<ARGS>::operator++() 
{    // Prefix ++
    ++ii;
    while (!ii) {
        ++k;
        if (k == M.outerSize()) break;    // Iteration is over
        ii.~InnerIterator();
        new (&ii) typename Eigen::SparseMatrix<ARGS>::InnerIterator(M,k);
    }
    return *this;
}


// -----------------------------------

template<class _Scalar, int _Options, class _StorageIndex>
EigenSparseMatrixIterator<ARGS> begin(Eigen::SparseMatrix<ARGS> const &M)
    { return EigenSparseMatrixIterator<ARGS>(M,0); }

template<class _Scalar, int _Options, class _StorageIndex>
EigenSparseMatrixIterator<ARGS> end(Eigen::SparseMatrix<ARGS> const &M)
    { return EigenSparseMatrixIterator<ARGS>(M,M.outerSize()); }

#undef ARGS
// -----------------------------------
template<class _Scalar, int _Options, class _StorageIndex>
#define ARGS _Scalar,_Options,_StorageIndex

blitz::Array<_Scalar,1> sum(
    Eigen::SparseMatrix<ARGS> const &M, int dimi, char invert='+');

/** Sum the rows or columns of an Eigen SparseMatrix.
@param dimi 0: sum rows, 1: sum columns */
template<class _Scalar, int _Options, class _StorageIndex>
blitz::Array<_Scalar,1> sum(
    Eigen::SparseMatrix<ARGS> const &M, int dimi, char invert='+')
{
    // Get our weight vector (in dense coordinate space)
    blitz::Array<_Scalar,1> ret;
    spcopy(
        accum::permute(accum::in_rank<2>(), {dimi},
        accum::blitz_new(ret)),
        M);

    if (invert == '-')
        for (int i=0; i<ret.extent(0); ++i) ret(i) = 1./ret(i);

    return ret;
}


template<class _Scalar, int _Options, class _StorageIndex>
Eigen::SparseMatrix<ARGS> sum_to_diagonal(
    Eigen::SparseMatrix<ARGS> const &M, int dimi, char invert='+');

template<class _Scalar, int _Options, class _StorageIndex>
Eigen::SparseMatrix<ARGS> sum_to_diagonal(
    Eigen::SparseMatrix<ARGS> const &M, int dimi, char invert='+')
{
    // Set up diagonal matrix w/ predictable indices in a TupleList
    _StorageIndex n = (dimi == 0 ? M.rows() : M.cols());
    TupleList<_StorageIndex, _Scalar, 2> tuples;
    tuples.set_shape({n,n});
    tuples.reserve(n);
    for (int i=0; i<n; ++i) tuples.add({i,i},0);

    for (auto ii(begin(M)); ii != end(M); ++ii) {
        auto ix(dimi == 0 ? ii->row() : ii->col());
        tuples[ix].value() += ii->value();
    }

    if (invert == '-') {
        for (size_t i=0; i<tuples.size(); ++i)
            tuples[i].value() = 1. / tuples[i].value();
    }

    return to_eigen(tuples);
}
#undef ARGS
// --------------------------------------------------------------
// --------------------------------------------------------------
template<class SparseIndexT, class _Scalar, int _Options, class _StorageIndex>
/** Create Eigen matrices with dense index space; while still
    assembling the dense index space.  Delay final creation of the
    matrix until all other matrices on the same index space are
    complete. */
class MakeDenseEigen {
public:
    typedef SparseSet<SparseIndexT, _StorageIndex> SparseSetT;

    template<int RANK>
        using TupleListT = TupleList<_StorageIndex,_Scalar,RANK>;

    typedef accum::Sparsify<
        accum::Permute<
            accum::Ref<
                TupleListT<2>>,
            2>,
        SparseSetT, typename SparseSetT::sparse_type> AccumT;
    typedef Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> EigenSparseMatrixT;

    TupleListT<2> M;
    std::array<SparseSetT *,2> const dims;
    std::array<int,2> const permute;    // {1,0} for transpse=='T'
    AccumT accum;

    /** @param transpose Use '.' for regular, 'T' for transpose.  If transposing
        in the final TupleList output, the SparseSets to which indices
        are added are NOT transposed.
    */
    MakeDenseEigen(
        std::function<void (AccumT &)> const &fn,
        SparsifyTransform sparsify_transform,
        std::array<SparseSetT *,2> const &_dims,
        char transpose = '.')    // '.' or 'T
    : dims(_dims),
        permute((transpose == 'T') ? ibmisc::make_array(1,0) : ibmisc::make_array(0,1)),
        accum(
        accum::sparsify(sparsify_transform,
            accum::in_index_type<typename SparseSetT::sparse_type>(),
            dims,
        accum::permute(accum::in_rank<2>(), permute,
        accum::ref(M))))
    {
        fn(accum);
    }

    EigenSparseMatrixT to_eigen()
    {
        Eigen::SparseMatrix<typename AccumT::val_type,0,typename AccumT::index_type> M(
            accum.dim(permute[0]).extent(),
            accum.dim(permute[1]).extent());

        // This segfaults if indices are out of bounds
        M.setFromTriplets(accum.base().begin(), accum.base().end());
        return M;
    }
};

// ------------------------------------------


// ====================================================
template<class _Scalar, int _Options, class _StorageIndex>
void nc_write_eigen(
    netCDF::NcGroup *nc,
    Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> *A,
    std::string const &vname);

template<class _Scalar, int _Options, class _StorageIndex>
void nc_write_eigen(
    netCDF::NcGroup *nc,
    Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> *A,
    std::string const &vname)
{
    netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
    netCDF::NcVar vals_v = nc->getVar(vname + ".values");

    std::vector<size_t> startp = {0, 0};        // SIZE, RANK
    std::vector<size_t> countp = {1, 2};  // Write RANK elements at a time
    for (auto ii = begin(*A); ii != end(*A); ++ii, ++startp[0]) {
        auto index(ii->index());
        auto &val(ii->value());

        indices_v.putVar(startp, countp, &index[0]);
        vals_v.putVar(startp, countp, &val);
    }
}



template<class _Scalar, int _Options, class _StorageIndex>
void ncio_eigen(
    ibmisc::NcIO &ncio,
    Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> &A,
    std::string const &vname);

template<class _Scalar, int _Options, class _StorageIndex>
void ncio_eigen(
    ibmisc::NcIO &ncio,
    Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> &A,
    std::string const &vname)
{
    if (ncio.rw == 'r') (*ibmisc::ibmisc_error)(-1,
        "ncio_eigen() currently does not support reading.");

    std::vector<std::string> const dim_names({vname + ".size", vname + ".rank"});
    std::vector<netCDF::NcDim> dims;        // Dimensions in NetCDF
    std::vector<size_t> dim_sizes;          // Length of our two dimensions.


    // Count the number of elements in the sparse matrix.
    // NOTE: This can/does give a diffrent answer from A.nonZeros().
    //       But it is what we want for dimensioning netCDF arrays.
    long count=0;
    for (auto ii(begin(A)); ii != end(A); ++ii) ++count;
    dims = ibmisc::get_or_add_dims(ncio, dim_names, {count, 2});

    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    std::array<size_t,2> shape { A.rows(), A.cols() };
    ibmisc::get_or_put_att(info_v, 'w', "shape", "int64", shape);

    get_or_add_var(ncio, vname + ".indices", "int64", dims);
    get_or_add_var(ncio, vname + ".values", "double", {dims[0]});
    ncio += std::bind(&nc_write_eigen<_Scalar, _Options, _StorageIndex>, ncio.nc, &A, vname);

}
// ====================================================
/** Reference Eigen::Matrix to blitz::Array
The Eigen::Matrix Rvalue-reference is stored in memory allocated by TmpAlloc,
which must live at least as long as the returned blitz::Array. */
template<class val_type>
blitz::Array<val_type,2> to_blitz(
    Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> &&M,
    ibmisc::TmpAlloc &tmp)
{
   auto &M_tmp(tmp.move<Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic>>(std::move(M)));

   return blitz::Array<double,2>(
       M_tmp.data(),
       blitz::shape(M_tmp.cols(), M_tmp.rows()),
       blitz::neverDeleteData);
}

/** Reference Eigen::Matrix vector to blitz::Array
The Eigen::Matrix Rvalue-reference is stored in memory allocated by TmpAlloc,
which must live at least as long as the returned blitz::Array. */
template<class val_type>
blitz::Array<val_type,1> to_blitz(
    Eigen::Matrix<val_type, Eigen::Dynamic, 1> &&M,
    ibmisc::TmpAlloc &tmp)
{
   auto &M_tmp(tmp.move<Eigen::Matrix<val_type, Eigen::Dynamic, 1>>(std::move(M)));

   return blitz::Array<double,1>(
       M_tmp.data(),
       blitz::shape(M_tmp.rows()),
       blitz::neverDeleteData);
}

/** Reference Eigen::Matrix vector to blitz::Array
The Eigen::Matrix Rvalue-reference is stored in memory allocated by TmpAlloc,
which must live at least as long as the returned blitz::Array. */
template<class val_type>
blitz::Array<val_type,1> to_blitz(
    Eigen::Matrix<val_type, 1, Eigen::Dynamic> &&M,
    ibmisc::TmpAlloc &tmp)
{
   auto &M_tmp(tmp.move<Eigen::Matrix<val_type, 1, Eigen::Dynamic>>(std::move(M)));

   return blitz::Array<double,1>(
       M_tmp.data(),
       blitz::shape(M_tmp.cols()),
       blitz::neverDeleteData);
}

// ==============================================================
// --------------------------------------------------------
template<class TupleT>
void consolidate(std::vector<TupleT> &A, bool zero_nan)
{

    // Sort the Tuples
    std::stable_sort(A.begin(), A.end());

    // Eliminate duplicates
    enum {NORUN, INRUN} state = NORUN;

    size_t dst=0;        // Where we're writing
    size_t start = 0;    // Start of run
    for (size_t i=0; i<A.size(); ++i) {
        switch(state) {
            case NORUN : {
                if (!isnone(A[i])) {
                    start = i;
                    state = INRUN;
                    A[dst] = A[i];
                }
            } break;
            case INRUN : {
                if (A[i].index() == A[start].index()) {
                    // Skip 0 and NaN
                    if (isnone(A[i].index(), zero_nan)) continue;

                    // Continue a run
                    A[dst].value += A[i].value;
                } else {
                    // Finish this run
                    ++dst;
                    state = NORUN;
                }
            } break;
        }
    }
}


/** @} */

}   // namespace

