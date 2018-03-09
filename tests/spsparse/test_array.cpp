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

// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <cstddef>
#include <functional>
#include <gtest/gtest.h>
#include <ibmisc/array.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/memory.hpp>
#include <spsparse/eigen.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/SparseSet.hpp>
#include <iostream>
#include <everytrace.h>

using namespace spsparse;
using namespace ibmisc;
using namespace std::placeholders;

// The fixture for testing class Foo.
class SpSparseTest : public ::testing::Test {
protected:
    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    SpSparseTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~SpSparseTest()
    {
        for (auto ii(tmpfiles.begin()); ii != tmpfiles.end(); ++ii) {
            ::remove(ii->c_str());
        }
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() {}

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown() {}

//    // The mock bar library shaed by all tests
//    MockBar m_bar;
};


TEST_F(SpSparseTest, TupleList) {
    // Make a simple TupleList
    TupleList<int, double, 1> arr1({4});
    arr1.add({1}, 2.);
    arr1.add({3}, 6.);
    EXPECT_EQ(2, arr1.size());
    {
        auto ii(arr1.begin());
        EXPECT_EQ(1, ii->index(0));
        EXPECT_EQ(2., ii->value());
        ++ii;
        EXPECT_EQ(3, ii->index(0));
    }

    // Test bounds checking
    try {
        arr1.add({17}, 4.);
        FAIL() << "Expected spsparse::Exception";
    } catch(ibmisc::Exception const &err) {
    } catch(...) {
        FAIL() << "Expected spsparse::Exception";
    }

    // Test Move Constructor
    TupleList<int, double, 1> arr2(std::move(arr1));
    EXPECT_EQ(0, arr1.size());
    EXPECT_EQ(2, arr2.size());
    {
        auto ii(arr2.begin());
        EXPECT_EQ(1, ii->index(0));
        EXPECT_EQ(2., ii->value());
        ++ii;
        EXPECT_EQ(3, ii->index(0));
    }
}

TEST_F(SpSparseTest, sort_tuple_list) {
    // Make a simple TupleList
    TupleList<int, double, 1> arr1({4});
    arr1.add({1}, 2.);
    arr1.add({3}, 6.);

    std::sort(arr1.tuples.begin(), arr1.tuples.end());
}

TEST_F(SpSparseTest, dense)
{
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr2({20,10});

    arr2.add({1,0}, 15.);
    arr2.add({1,3}, 17.);
    arr2.add({2,4}, 17.);
    arr2.add({6,4}, 10.);

    blitz::Array<double, 2> dense(spsparse::to_blitz(arr2));
    int i,j;
    double sum=0;
    for (int i=0; i<20; ++i) {
    for (int j=0; j<10; ++j) {
        sum += dense(i,j);
    }}
    EXPECT_EQ(sum, 59.);

    for (auto ii(arr2.begin()); ii != arr2.end(); ++ii)
        EXPECT_EQ(dense(ii->index(0), ii->index(1)), ii->value());

}


TEST_F(SpSparseTest, dense_to_blitz)
{
    typedef TupleList<int, double, 2> TupleListT;

    blitz::Array<double,2> dense1(4,5);
    dense1 = 0;
    dense1(2,3) = 5.0;
    dense1(2,4) = 6.0;
    dense1(0,1) = 7.0;

    TupleListT sparse1({4,5});
    TupleListT sparse2({4,5});
    spcopy(accum::tee(accum::ref(sparse1), accum::ref(sparse2)), dense1);    // Test multiple items for accum::refs()
    auto dense2(spsparse::to_blitz(sparse1));
    auto dense3(spsparse::to_blitz(sparse2));

    for (int i=0; i<dense1.extent(0); ++i) {
    for (int j=0; j<dense1.extent(0); ++j) {
        double const d1(dense1(i,j));
        double const d2(dense2(i,j));
        double const d3(dense2(i,j));

        EXPECT_EQ(d1, d2);
        EXPECT_EQ(d1, d3);
    }}
}

TEST_F(SpSparseTest, sparse_set_ncio)
{
    // If the constructor and destructor are not enough for setting up
    std::string fname("__sparse_set_test.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());

    SparseSet<long,int> dim0;
    dim0.add_dense(17);
    dim0.add_dense(3);
    dim0.add_dense(5000);

    {ibmisc::NcIO ncio(fname, NcFile::replace);
        dim0.ncio(ncio, "dim0");
    }

    SparseSet<long,int> dim1;
    {ibmisc::NcIO ncio(fname, NcFile::read);
        dim1.ncio(ncio, "dim0");
    }

    EXPECT_EQ(dim0, dim1);
}

TEST_F(SpSparseTest, sparse_set)
{
    // Construct a SparseMatrix
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr2({20,10});

    arr2.add({6,4}, 10.);
    arr2.add({1,0}, 15.);
    arr2.add({2,4}, 17.);
    arr2.add({1,3}, 17.);

    // Make a SparseSet for each of its dimensions
    SparseSet<int,int> dim0;
    dim0.set_sparse_extent(arr2.shape()[0]);
    dim0.add_sorted(
        dim_index_iter(arr2.begin(), 0),
        dim_index_iter(arr2.end(), 0));

    EXPECT_EQ(3, dim0.dense_extent());
    EXPECT_EQ(1, dim0.to_sparse(0));
    EXPECT_EQ(2, dim0.to_sparse(1));
    EXPECT_EQ(6, dim0.to_sparse(2));
    EXPECT_EQ(0, dim0.to_dense(1));
    EXPECT_EQ(1, dim0.to_dense(2));
    EXPECT_EQ(2, dim0.to_dense(6));

    // Test to_dense
    TupleListT arr2d;
    std::array<SparseSet<int,int> *, 2> dims({&dim0, nullptr});
    spcopy(
        accum::sparsify(
            accum::sparsify_normalize_transforms({SparsifyTransform::TO_DENSE}, dims),
            accum::in_index_type<int>(),
            dims,
        accum::ref(arr2d)),
        arr2);

    {
    auto ii(arr2d.begin());
    EXPECT_EQ(2, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    ++ii;
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    ++ii;
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(3, ii->index(1));
    ++ii;
    EXPECT_EQ(ii, arr2d.end());
    }

    // Test to_sparse
    {
    TupleListT arr3;
    spcopy(
        accum::to_sparse(
            make_array(&dim0, nullptr),
        accum::ref(arr3)),
        arr2d);

    auto ii(arr3.begin());
    EXPECT_EQ(6, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    ++ii;
    EXPECT_EQ(2, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(3, ii->index(1));
    ++ii;
    EXPECT_EQ(ii, arr3.end());
    }

}

TEST_F(SpSparseTest, partial_sparsify)
{
    // Construct a SparseMatrix
    TupleList<long, double, 2> arr({30,10});

    arr.add({6,4}, 10.);
    arr.add({17,0}, 15.);
    arr.add({22,4}, 17.);
    arr.add({17,3}, 17.);

    SparseSet<long,int> dim0;
    SparseSet<long,int> dim1;
    TupleList<int, double, 2> arr2;
    spcopy(
        accum::add_dense(ibmisc::make_array(&dim0, nullptr),
        accum::ref(arr2)),
        arr, false);

    auto ii(arr2.begin());
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    EXPECT_EQ(10., ii->value());
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    EXPECT_EQ(15., ii->value());
    ++ii;
    EXPECT_EQ(2, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    EXPECT_EQ(17., ii->value());
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(3, ii->index(1));
    EXPECT_EQ(17., ii->value());
    ++ii;

    EXPECT_TRUE(ii == arr2.end());
}

TEST_F(SpSparseTest, invert_accum)
{
    // Construct a SparseMatrix
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr({20,10});

    arr.add({6,4}, 8.);    // Inverses of powers of 2 will be exact FP arithmetic
    arr.add({1,0}, 16.);

    TupleListT arr2;
    spcopy(
        accum::invert('-',
        accum::ref(arr2)),
        arr);

    auto ii(arr2.begin());
    EXPECT_EQ(6, ii->index(0));
    EXPECT_EQ(4, ii->index(1));
    EXPECT_EQ(1./8., ii->value());
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    EXPECT_EQ(1./16., ii->value());
}


TEST_F(SpSparseTest, transpose_accum)
{
    // Construct a SparseMatrix
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr({20,10});

    arr.add({6,4}, 8.);    // Inverses of powers of 2 will be exact FP arithmetic
    arr.add({1,0}, 16.);

    TupleListT arr2;
    spcopy(accum::transpose(accum::ref(arr2)), arr);

    EXPECT_EQ(arr2.shape()[0], 10);
    EXPECT_EQ(arr2.shape()[1], 20);
    auto ii(arr2.begin());

    EXPECT_EQ(4, ii->index(0));
    EXPECT_EQ(6, ii->index(1));
    EXPECT_EQ(8., ii->value());
    ++ii;
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(1, ii->index(1));
    EXPECT_EQ(16., ii->value());

}
// -----------------------------------------------
TEST_F(SpSparseTest, sparse_set_accum)
{
    SparseSet<long,int> dimA, dimB;

    accum::SparseSetAccum<SparseSet<long,int>, double, 3> acc(
        {&dimA, nullptr, &dimB});

    acc.add({6L,7L, 2L}, 8.);
    acc.add({1L,7L, 1L}, 16.);
    acc.add({6L,8L, 0L}, 17.);

    EXPECT_EQ(2, dimA.dense_extent());
    EXPECT_EQ(true, dimA.in_sparse(6L));
    EXPECT_EQ(true, dimA.in_sparse(1L));
    EXPECT_EQ(6L, dimA.to_sparse(0));
    EXPECT_EQ(1L, dimA.to_sparse(1));

    EXPECT_EQ(3, dimB.dense_extent());
    EXPECT_EQ(true, dimB.in_sparse(2L));
    EXPECT_EQ(true, dimB.in_sparse(1L));
    EXPECT_EQ(true, dimB.in_sparse(0L));
    EXPECT_EQ(2L, dimB.to_sparse(0));
    EXPECT_EQ(1L, dimB.to_sparse(1));
    EXPECT_EQ(0L, dimB.to_sparse(2));
}
// -----------------------------------------------
typedef MakeDenseEigen<long,double,0,int> MakeDenseEigenT;

void sample_makedense(MakeDenseEigenT::AccumT &&accum)
{
    accum.add({6L,2L}, 8.);
    accum.add({1L,0L}, 16.);
    accum.add({6L,0L}, 17.);
}
TEST_F(SpSparseTest, make_dense_eigen)
{
    std::array<MakeDenseEigenT::SparseSetT, 2> dims;
    MakeDenseEigenT M_m(
        std::bind(&sample_makedense, _1),
        {SparsifyTransform::ADD_DENSE},
        {&dims[0], &dims[1]}, '.');

    auto ii(M_m.M.begin());
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    EXPECT_EQ(8., ii->value());
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(1, ii->index(1));
    EXPECT_EQ(16., ii->value());
    ++ii;
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(1, ii->index(1));
    EXPECT_EQ(17., ii->value());

    EXPECT_EQ(-1, M_m.M.shape(0));
    EXPECT_EQ(-1, M_m.M.shape(1));

    auto M(M_m.to_eigen());

    EXPECT_EQ(2, M.rows());
    EXPECT_EQ(2, M.cols());

    // Test simple iterator for Eigen::SparseMatrix
    {
    auto ii(begin(M));
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    EXPECT_EQ(8., ii->value());
    ++ii;
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(1, ii->index(1));
    EXPECT_EQ(17., ii->value());
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(1, ii->index(1));
    EXPECT_EQ(16., ii->value());
    ++ii;
    EXPECT_TRUE(ii == end(M));
    }
}

#if 0
TEST_F(SpSparseTest, diag_matrix)
{
    blitz::Array<double,1> diag(5);
    diag(0) = 4.;
    diag(1) = 8.;
    diag(2) = 0.;
    diag(3) = 16.;
    diag(4) = 0.;

    Eigen::DiagonalMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Md(map_eigen_diagonal(diag));
    Eigen::SparseMatrix<double> M(Md);

    EXPECT_EQ(5, M.rows());
    EXPECT_EQ(5, M.cols());


    auto ii(begin(M));
    EXPECT_EQ(0, ii->index(0));
    EXPECT_EQ(0, ii->index(1));
    EXPECT_EQ(4., ii->value());
    ++ii;
    EXPECT_EQ(1, ii->index(0));
    EXPECT_EQ(1, ii->index(1));
    EXPECT_EQ(8., ii->value());
    ++ii;
    EXPECT_EQ(3, ii->index(0));
    EXPECT_EQ(3, ii->index(1));
    EXPECT_EQ(16., ii->value());
    ++ii;
    EXPECT_TRUE(ii == end(M));
}
#endif

/** Not strictly a test of Eigen stuff.  But this is used in IceBin/RegridMatrices.cpp */
TEST_F(SpSparseTest, diag_matrix_inv)
{
    blitz::Array<double,1> diag(4);
    diag(0) = 4.;
    diag(1) = 8.;
    diag(2) = 0.;
    diag(3) = 16.;

    blitz::Array<double,1> diag_inv(1. / diag);

    EXPECT_EQ(1./4., diag_inv(0));
    EXPECT_EQ(1./8., diag_inv(1));
    EXPECT_EQ(1./16., diag_inv(3));
}


TEST_F(SpSparseTest, sum)
{
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr({3,4});

    // Note there is NOTHING in column 0 of the matrix.
    // This tests correctness of the begin() and end() iterators
    arr.add({0,2}, 1.);
    arr.add({0,1}, 2.);
    arr.add({0,3}, 4.);

    arr.add({2,3}, 8.);
    arr.add({2,1}, 16.);

    Eigen::SparseMatrix<double,0,int> M(arr.shape(0), arr.shape(1));
    M.setFromTriplets(arr.begin(), arr.end());

    // Test correctness of iterator
    int n=0;
    for (auto ii(begin(M)); ii != end(M); ++ii) {
        // printf("M(%d, %d) = %f\n", ii->row(), ii->col(), ii->value());
        ++n;
    }
    EXPECT_EQ(5, n);


    {auto arrb(sum(M,0,'+'));
        EXPECT_EQ(3, arrb.extent(0));
        EXPECT_EQ(7., arrb(0));
        EXPECT_EQ(0., arrb(1));
        EXPECT_EQ(24., arrb(2));
    }

    {auto arrb(sum(M,0,'-'));
        EXPECT_EQ(3, arrb.extent(0));
        EXPECT_DOUBLE_EQ(1./7., arrb(0));
//        EXPECT_EQ(0., arrb(1));
        EXPECT_DOUBLE_EQ(1./24., arrb(2));
    }

    {auto arrb(sum(M,1,'+'));
        EXPECT_EQ(4, arrb.extent(0));
        EXPECT_EQ(0., arrb(0));
        EXPECT_EQ(18., arrb(1));
        EXPECT_EQ(1., arrb(2));
        EXPECT_EQ(12., arrb(3));
    }
}

// ------------------------------------------------------------------
void *sample_eigen_data;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> sample_eigen_to_blitz()
{
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M(2,3);
    int count=0;

    for (int j=0; j<3; ++j) {
    for (int i=0; i<2; ++i) {
        M(i,j) = count++;
    }}

    sample_eigen_data = M.data();
    return M;
}

TEST_F(SpSparseTest, eigen_to_blitz)
{
    TmpAlloc tmp;
    blitz::Array<double,2> M(to_blitz(sample_eigen_to_blitz(), tmp));

    // Test that the object created in sample_eigen_to_blitz()
    // is moved, not copied, to the TmpAlloc
    EXPECT_TRUE(sample_eigen_data == M.data());
    EXPECT_EQ(3, M.extent(0));
    EXPECT_EQ(2, M.extent(1));

    // Test that it has the value we expect (transposing indices
    // for Eigen vs. Blitz)
    int count=0;
    for (int j=0; j<M.extent(0); ++j) {
    for (int i=0; i<M.extent(1); ++i) {    // stride=1
        EXPECT_EQ((double)count, M(j,i));
        ++count;
    }}
}
// ------------------------------------------------------------------
Eigen::Matrix<double,Eigen::Dynamic,1> sample_eigen_to_blitz_1d()
{
    Eigen::Matrix<double,Eigen::Dynamic,1> M(3);
    int count=0;

    for (int i=0; i<3; ++i) {
        M[i] = count++;
    }

    sample_eigen_data = M.data();
    return M;
}

TEST_F(SpSparseTest, eigen_to_blitz_1d)
{
    TmpAlloc tmp;
    blitz::Array<double,1> M(to_blitz(sample_eigen_to_blitz_1d(), tmp));

    // Test that the object created in sample_eigen_to_blitz()
    // is moved, not copied, to the TmpAlloc
    EXPECT_TRUE(sample_eigen_data == M.data());
    EXPECT_EQ(3, M.extent(0));

    // Test that it has the value we expect (transposing indices
    // for Eigen vs. Blitz)
    int count=0;
    for (int i=0; i<M.extent(1); ++i) {
        EXPECT_EQ((double)count, M(i));
        ++count;
    }
}
// ------------------------------------------------------------------
TEST_F(SpSparseTest, reshape)
{
    blitz::Array<double,2> a2(2,3);

    auto a3(reshape<double,2,2>(a2, {3,2}));
    EXPECT_TRUE(a2.data() == a3.data());
    EXPECT_EQ(2, a3.rank());
    EXPECT_EQ(3, a3.extent(0));
    EXPECT_EQ(2, a3.extent(1));
    EXPECT_EQ(1, a3.stride()[1]);
    EXPECT_EQ(2, a3.stride()[0]);

    auto a1(reshape<double,2,1>(a2, {-1}));
    EXPECT_TRUE(a1.data() == a3.data());
    EXPECT_EQ(1, a1.rank());
    EXPECT_EQ(6, a1.extent(0));
    EXPECT_EQ(1, a1.stride()[0]);

    EXPECT_ANY_THROW(({
        reshape<double,2,2>(a2, {-1,-1});
    }));
    EXPECT_ANY_THROW(({
        reshape<double,2,2>(a2, {5,3});
    }));
    EXPECT_ANY_THROW(({
        reshape<double,2,2>(a2, {5,-1});
    }));
    EXPECT_ANY_THROW(({
        reshape<double,2,2>(a2, {-1,5});
    }));

}
// ------------------------------------------------------------------
TEST_F(SpSparseTest, indices)
{
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr({3,4});

    arr.add({0,0}, 1.);
    arr.add({2,3}, 8.);
    arr.add({2,3}, 9.);

    std::vector<std::array<int,2>> idx;
    spcopy(
        accum::Indices<int,double,2>(idx),
        arr);

    EXPECT_EQ(3, idx.size());
    EXPECT_TRUE((idx[0] == std::array<int,2>{0,0}));
    EXPECT_TRUE((idx[1] == std::array<int,2>{2,3}));
    EXPECT_TRUE((idx[2] == std::array<int,2>{2,3}));

}
// ------------------------------------------------------------------


int main(int argc, char **argv) {
    // For test, configure Everytrace to silently throw exceptions (which we can catch)
    everytrace_init();
    everytrace_exit = &everytrace_exit_silent_exception;

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
