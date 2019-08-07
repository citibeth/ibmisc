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
#include <ibmisc/error.hpp>
#include <ibmisc/iter.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/blitz.hpp>
#include <spsparse/SparseSet.hpp>
#include <spsparse/tuplelist.hpp>

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
    bool set_shape)
{
    if (set_shape) ret.set_shape(std::array<long,2>{M.rows(), M.cols()});

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

// --------------------------------------------------------------

// --------------------------------------------------------
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

    // Run this every time ii is changed.
    void validate_ii()
    {
        while (k != M.outerSize() && !ii) {
            ++k;
            ii.~InnerIterator();
            new (&ii) typename Eigen::SparseMatrix<ARGS>::InnerIterator(M,k);
        }
    }
public:
    typedef _StorageIndex index_type;
    typedef _Scalar val_type;
    static const int RANK = 2;

    Eigen::SparseMatrix<ARGS> const &M;
    int k;
    typename Eigen::SparseMatrix<ARGS>::InnerIterator ii;
    
    EigenSparseMatrixIterator(Eigen::SparseMatrix<ARGS> const &_M, int _k)
        : M(_M), k(_k),
        ii(typename Eigen::SparseMatrix<ARGS>::InnerIterator(M,k)) {
            validate_ii();
        }

    // ---------------------------------------------------------
    // http://www.cplusplus.com/reference/iterator/

    IteratorT &operator++();

    bool operator==(IteratorT const &other) const
    {
        if (k != other->k) return false;
        return (k == M.outerSize()) || (ii == other->ii);
    }

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
    validate_ii();
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
@param dimi The dimension to keep.
   0: sum columns, leaving rows (a column vector)
   1: sum rows, leaving columns (a row vector) */
template<class _Scalar, int _Options, class _StorageIndex>
blitz::Array<_Scalar,1> sum(
    Eigen::SparseMatrix<ARGS> const &M, int dimi, char invert)
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
std::array<blitz::Array<_Scalar,1>,2> sums(
    Eigen::SparseMatrix<ARGS> const &M, char invert='+');

/** Sum the rows or columns of an Eigen SparseMatrix.
@param dimi 0: sum rows, 1: sum columns */
template<class _Scalar, int _Options, class _StorageIndex>
std::array<blitz::Array<_Scalar,1>,2> sums(
    Eigen::SparseMatrix<ARGS> const &M, char invert)
{
    // Get our weight vector (in dense coordinate space)
    std::array<blitz::Array<_Scalar,1>,2> rets;
    rets[0].reference(blitz::Array<_Scalar,1>(M.rows()));
    rets[1].reference(blitz::Array<_Scalar,1>(M.cols()));

    rets[0] = 0;
    rets[1] = 0;

    // See here for note on iterating through Eigen::SparseMatrix
    // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1104
    for (int k=0; k<M.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<ARGS>::InnerIterator ii(M,k); ii; ++ii) {
        rets[0](ii.row()) += ii.value();
        rets[1](ii.col()) += ii.value();
    }}

    if (invert == '-') {
        for (int j=0; j<2; ++j) {
        for (int i=0; i<rets[j].extent(0); ++i) {
            rets[j](i) = 1./rets[j](i);
        }}
    }

    return rets;
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

    typedef
        accum::Sparsify<
        accum::Permute<
            accum::IncludeZero<
            accum::Ref<
                TupleListT<2>>>,
            2>,
        SparseSetT, typename SparseSetT::sparse_type> AccumT;
    typedef Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> EigenSparseMatrixT;

    TupleListT<2> M;
    std::array<SparsifyTransform,2> const normalized_transform;
    std::array<SparseSetT *,2> const dims;
    std::array<int,2> const permute;    // {1,0} for transpse=='T'
    bool const include_zero;    // Include zero elements in the matrix?  (Or filter out...)

    /** Constructs and returns the accumulator appropriate for this
        matrix generator. */
    AccumT accum()
    {
        return
            accum::sparsify(normalized_transform,
                accum::in_index_type<typename SparseSetT::sparse_type>(),
                dims,
            accum::permute(accum::in_rank<2>(), permute,
            accum::include_zero(include_zero,
            accum::ref(M))));
    }

    /** @param transpose Use '.' for regular, 'T' for transpose.  If transposing
        in the final TupleList output, the SparseSets to which indices
        are added are NOT transposed.
    */
    MakeDenseEigen(
        std::vector<SparsifyTransform> const &_sparsify_transform,
        std::array<SparseSetT *,2> const &_dims,
        char transpose = '.',    // '.' or 'T
        bool _include_zero=true);

    MakeDenseEigen(
        std::function<void (AccumT &&)> const &fn,
        std::vector<SparsifyTransform> const &_sparsify_transform,
        std::array<SparseSetT *,2> const &_dims,
        char transpose = '.',    // '.' or 'T
        bool _include_zero=true)
    : MakeDenseEigen(_sparsify_transform, _dims, transpose, _include_zero)
    {
        fn(accum());
    }


    long extent(int i)
    {
        switch(normalized_transform[i]) {
            case SparsifyTransform::ID :
                return -1;    // sparse_set == NULL
            case SparsifyTransform::ADD_DENSE :
            case SparsifyTransform::TO_DENSE_IGNORE_MISSING :
            case SparsifyTransform::TO_DENSE :
            case SparsifyTransform::KEEP_DENSE :
                return dims[i]->dense_extent();
            case SparsifyTransform::TO_SPARSE :
            case SparsifyTransform::KEEP_SPARSE :
                return dims[i]->sparse_extent();

        }
    }

    EigenSparseMatrixT to_eigen();
};


template<class SparseIndexT, class _Scalar, int _Options, class _StorageIndex>
MakeDenseEigen<SparseIndexT,_Scalar,_Options,_StorageIndex>::MakeDenseEigen(
    std::vector<SparsifyTransform> const &_sparsify_transform,
    std::array<SparseSetT *,2> const &_dims,
    char transpose,    // '.' or 'T
    bool _include_zero)
: normalized_transform(spsparse::accum::sparsify_normalize_transforms(_sparsify_transform, _dims)),
    dims(_dims),
    include_zero(_include_zero),
    permute((transpose == 'T') ? ibmisc::make_array(1,0) : ibmisc::make_array(0,1))
{}

#define ARGS SparseIndexT,_Scalar,_Options,_StorageIndex
template<class SparseIndexT, class _Scalar, int _Options, class _StorageIndex>
typename MakeDenseEigen<ARGS>::EigenSparseMatrixT MakeDenseEigen<ARGS>::to_eigen()
{
    for (int i=0; i<2; ++i) {
        if (extent(i) < 0) {
            (*ibmisc::ibmisc_error)(-1, "MakeDenseEigen requires dimensions to have a computed extent.  It is not compatible with SparsifyTransform::ID");
        }
    }

    Eigen::SparseMatrix<typename AccumT::val_type,0,typename AccumT::index_type> Matrix(
        extent(permute[0]),
        extent(permute[1]));

    // This segfaults if indices are out of bounds
    Matrix.setFromTriplets(M.begin(), M.end());
    return Matrix;
}
#undef ARGS



// ------------------------------------------


// ====================================================
template<class _Scalar, int _Options, class _StorageIndex>
void nc_rw_eigen(
    netCDF::NcGroup *nc,
    char rw,
    Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> *A,
    std::string const &vname);

template<class _Scalar, int _Options, class _StorageIndex>
static void nc_rw_eigen(
    netCDF::NcGroup *nc,
    char rw,
    Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> *A,
    std::string const &vname)
{
    netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
    netCDF::NcVar vals_v = nc->getVar(vname + ".values");

    if (rw == 'w') {
        int const N = A->nonZeros();

        // Create in-memory data structure amenable to writing to disk quickly
        {std::vector<int> indices;
            indices.reserve(N*2);
            for (auto ii = begin(*A); ii != end(*A); ++ii) {
                indices.push_back(ii->row());
                indices.push_back(ii->col());
            }
            indices_v.putVar(&indices[0]);
        }

        // Write it out!
        {std::vector<double> vals;
            vals.reserve(N);
            for (auto ii = begin(*A); ii != end(*A); ++ii) {
                vals.push_back(ii->value());
            }
            vals_v.putVar(&vals[0]);    // Write to entire NetCDF variable directly from RAM
        }
    } else {    // rw == 'r'
        int const N = vals_v.getDim(0).getSize();

        // Create in-memory data structure amenable to writing to disk quickly
        std::vector<int> indices;
        std::vector<double> vals;

        // Read it
        indices.resize(N*2);
        indices_v.getVar(&indices[0]);
        vals.resize(N);
        vals_v.getVar(&vals[0]);    // Write to entire NetCDF variable directly from RAM

        // Read the shape
        ibmisc::NcIO ncio(nc, rw);    // dummy
        auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
        std::vector<size_t> shape;
        ibmisc::get_or_put_att(info_v, rw, "shape", "int64", shape);

        // Convert to TupleList
        // TODO: If we are more clever with iterators, we don't have to copy to
        //       a TupleList first.
        TupleList<_StorageIndex,_Scalar,2> tuples;
        for (size_t i=0; i<vals.size(); ++i) {
            tuples.add({indices[i*2], indices[i*2+1]}, vals[i]);
        }
        indices.clear();
        vals.clear();

        A->setZero();
        A->resize(shape[0], shape[1]);
        A->reserve(tuples.tuples.size());

        // Convert to Eigen
        A->setFromTriplets(tuples.begin(), tuples.end());
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
    std::vector<std::string> const dim_names({vname + ".nnz", vname + ".rank"});
    std::vector<netCDF::NcDim> dims;        // Dimensions in NetCDF

    // Count the number of elements in the sparse matrix.
    // NOTE: This can/does give a diffrent answer from A.nonZeros().
    //       But it is what we want for dimensioning netCDF arrays.
    if (ncio.rw == 'w') {
        long const count = (ncio.rw == 'w' ? A.nonZeros() : 0);
        dims = ibmisc::get_or_add_dims(ncio, dim_names, {count, 2});

        auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
        std::array<size_t,2> shape { A.rows(), A.cols() };
        ibmisc::get_or_put_att(info_v, 'w', "shape", "int64", shape);
    } else {    // Read 
        dims = ibmisc::get_or_add_dims(ncio, dim_names, {0,0});    // length ignored on read
    }

    get_or_add_var(ncio, vname + ".indices", "int", dims);
    get_or_add_var(ncio, vname + ".values", "double", {dims[0]});
    ncio += std::bind(&nc_rw_eigen<_Scalar, _Options, _StorageIndex>, ncio.nc, ncio.rw, &A, vname);

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
   auto &M_tmp(tmp.take<Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic>>(std::move(M)));

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
   auto &M_tmp(tmp.take<Eigen::Matrix<val_type, Eigen::Dynamic, 1>>(std::move(M)));

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
   auto &M_tmp(tmp.take<Eigen::Matrix<val_type, 1, Eigen::Dynamic>>(std::move(M)));

   return blitz::Array<double,1>(
       M_tmp.data(),
       blitz::shape(M_tmp.cols()),
       blitz::neverDeleteData);
}

// ==============================================================
// ------------------------------------------------------------------
template<class val_type>
blitz::Array<val_type,2> to_blitz(
    Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> &M)
{
   return blitz::Array<double,2>(
       M.data(),
       blitz::shape(M.cols(), M.rows()),
       blitz::neverDeleteData);
}

template<class val_type>
blitz::Array<val_type,2> const to_blitz(
    Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> const &M)
{
   return blitz::Array<double,2>(
       const_cast<val_type *>(M.data()),
       blitz::shape(M.cols(), M.rows()),
       blitz::neverDeleteData);
}
// ------------------------------------------------------------------
template<class val_type>
blitz::Array<val_type,1> to_blitz(
    Eigen::Matrix<val_type, Eigen::Dynamic, 1> &M)
{
   return blitz::Array<double,1>(
       M.data(),
       blitz::shape(M.rows()),
       blitz::neverDeleteData);
}

template<class val_type>
blitz::Array<val_type,1> const to_blitz(
    Eigen::Matrix<val_type, Eigen::Dynamic, 1> const &M)
{
   return blitz::Array<double,1>(
       const_cast<val_type *>(M.data()),
       blitz::shape(M.rows()),
       blitz::neverDeleteData);
}
// ------------------------------------------------------------------
template<class val_type>
blitz::Array<val_type,1> to_blitz(
    Eigen::Matrix<val_type, 1, Eigen::Dynamic> &M)
{
   return blitz::Array<double,1>(
       M.data(),
       blitz::shape(M.cols()),
       blitz::neverDeleteData);
}
template<class val_type>
blitz::Array<val_type,1> const to_blitz(
    Eigen::Matrix<val_type, 1, Eigen::Dynamic> &M)
{
   return blitz::Array<double,1>(
       const_cast<val_type *>(M.data()),
       blitz::shape(M.cols()),
       blitz::neverDeleteData);
}

// ==============================================================

template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>
    map_eigen_matrix(blitz::Array<Scalar,2> &A_b);

/** View 2-D blitz::Array as an Eigen::Matrix */
template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>
    map_eigen_matrix(blitz::Array<Scalar,2> &A_b)
{
    if (is_row_major(A_b)) {
        return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(
            A_b.data(), A_b.extent(1), A_b.extent(0));
    } else if (is_column_major(A_b)) {
        return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(
            A_b.data(), A_b.extent(0), A_b.extent(1));
    } else (*ibmisc::ibmisc_error)(-1, "Matrix must be row-major or column-major");
}

/** View 1-D blitz::Array as an Eigen column vector */
template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>
    map_eigen_colvector(blitz::Array<Scalar,1> &A_b)
{
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(
        A_b.data(), A_b.extent(0));
}
template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> const
    map_eigen_colvector(blitz::Array<Scalar,1> const &A_b)
{
    return map_eigen_colvector(
        const_cast<blitz::Array<Scalar,1> &>(A_b));
}



/** View 1-D blitz::Array as an Eigen row vector */
template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,1,Eigen::Dynamic>>
    map_eigen_rowvector(blitz::Array<Scalar,1> &A_b)
{
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(
        A_b.data(), A_b.extent(0));
}
template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> const
    map_eigen_rowvector(blitz::Array<Scalar,1> const &A_b)
{
    return map_eigen_rowvector(
        const_cast<blitz::Array<Scalar,1> &>(A_b));
}




/** View 1-D blitz::Array as a diagonal Eigen matrix */
template<class Scalar>
const Eigen::DiagonalWrapper<const Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>>
    map_eigen_diagonal(blitz::Array<Scalar,1> &A_b)
{
    return map_eigen_colvector(A_b).asDiagonal();
}
template<class Scalar>
Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> const
    map_eigen_diagonal(blitz::Array<Scalar,1> const &A_b)
{
    return map_eigen_diagonal(
        const_cast<blitz::Array<Scalar,1> &>(A_b));
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

