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

#include <memory>
#include <ibmisc/cython.hpp>
#include <spsparse/eigen.hpp>
#include <spsparse/SparseSet.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <ibmisc/enum.hpp>

#include "examples.hpp"

using namespace spsparse;
using namespace ibmisc;

namespace ibmisc {
namespace cython {

// -----------------------------------------
// Types that will be used throughout as template arguments
typedef long sparse_index_type;
typedef int dense_index_type;
typedef double val_type;

typedef spsparse::MakeDenseEigen<sparse_index_type, val_type, 0, dense_index_type> MakeDenseEigenT;
template<int RANK>
    using TupleListT = MakeDenseEigenT::TupleListT<RANK>;
template<int RANK>
    using DenseArrayT = blitz::Array<val_type,RANK>;
typedef MakeDenseEigenT::SparseSetT SparseSetT;
typedef MakeDenseEigenT::EigenSparseMatrixT EigenSparseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic> EigenRowVectorT;
typedef Eigen::DiagonalMatrix<val_type, Eigen::Dynamic> EigenDiagonalMatrixT;
// -----------------------------------------


static void example_double_blitz(blitz::Array<double,1> &a)
{
    for (int i=a.lbound(0); i<=a.ubound(0); ++i) {
        a(i) *= 2;
    }
}

/** Cython interface... converts np.array to blitz::Array */
void cyexample_double_blitz(PyObject *a)
{
    blitz::Array<double,1> ab(np_to_blitz<double,1>(a, "A", {-1}));
    example_double_blitz(ab);
}

PyObject *cyexample_sparse_matrix()
{
    spsparse::TupleList<long, double, 2> M({4,6});
    M.add({2,3}, 2.2);
    M.add({0,5}, 1.1);
    PyObject *tuple = spsparse_to_tuple(M);
    return tuple;
}

static void sample_makedense(MakeDenseEigenT::AccumT &&accum)
{

    accum.add({0,0},.5);
    accum.add({10,10},.5);
    accum.add({20,20},.5);


    //accum.add({30,30},.5);    // Nullspace
    accum.add({00,40},.5);
    accum.add({10,40},.5);
    accum.add({20,40},.5);
    //accum.add({30,40},.5);
}

std::unique_ptr<linear::Weighted> example_linear_weighted(
    std::string slinear_type)
{
    // ============ BvA1
    // Compute an Eigen-based overlap matrix
    std::array<MakeDenseEigenT::SparseSetT, 2> dims;

    dims[0].set_sparse_extent(40);
    dims[1].set_sparse_extent(50);

    // Add extra unneeded dimensions
    dims[0].add_dense(15);
    dims[1].add_dense(15);

    MakeDenseEigenT BvA_o(
        std::bind(&sample_makedense, std::placeholders::_1),
        {spsparse::SparsifyTransform::ADD_DENSE},
        {&dims[0], &dims[1]}, '.');

    // Add extra unneeded dimensions
    dims[0].add_dense(25);
    dims[1].add_dense(25);

 
    auto BvA(BvA_o.to_eigen());

    // Produce a scaled regridding matrix
    auto sBvA(sum(BvA, 0, '-'));

    // Scale and stick into a linear::Weighted_Eigen
    std::unique_ptr<linear::Weighted_Eigen> BvA1(new linear::Weighted_Eigen);
    BvA1->tmp.take(BvA1->dims[0], std::move(dims[0]));
    BvA1->tmp.take(BvA1->dims[1], std::move(dims[1]));
//    BvA1->dims[0] = &dims[0];
//    BvA1->dims[1] = &dims[1];
    BvA1->M.reset(new MakeDenseEigenT::EigenSparseMatrixT(map_eigen_diagonal(sBvA) * BvA));
    BvA1->wM.reference(sum(BvA,0,'+'));
    BvA1->Mw.reference(sum(BvA,1,'+'));

    // Make it not conservative
    BvA1->wM(2) *= .9;
    BvA1->conservative = false;

    auto linear_type(parse_enum<linear::LinearType>(slinear_type));

    if (linear_type == linear::LinearType::EIGEN) {
        auto ret(std::unique_ptr<linear::Weighted>(BvA1.release()));
        return ret;
    } else {
        auto ret(std::unique_ptr<linear::Weighted>(
            new linear::Weighted_Compressed(
                compress(*BvA1))));
        return ret;
    }
}


}}
