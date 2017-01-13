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
#include <gtest/gtest.h>
#include <ibmisc/array.hpp>
#include <spsparse/eigen.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/SparseSet.hpp>
#include <iostream>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

using namespace spsparse;
using namespace ibmisc;

// The fixture for testing class Foo.
class SpSparseTest : public ::testing::Test {
protected:

    // You can do set-up work for each test here.
    SpSparseTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~SpSparseTest() {}

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
        FAIL() << "Excpected spsparse::Exception";
    } catch(ibmisc::Exception const &err) {
    } catch(...) {
        FAIL() << "Excpected spsparse::Exception";
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
    spcopy(sparse1, dense1);
    auto dense2(spsparse::to_blitz(sparse1));

    for (int i=0; i<dense1.extent(0); ++i) {
    for (int j=0; j<dense1.extent(0); ++j) {
        double const d1(dense1(i,j));
        double const d2(dense2(i,j));

        EXPECT_EQ(d1, d2);
    }}

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
    spcopy(
        accum::sparse_transform(
            SparseTransform::TO_DENSE,
            ibmisc::make_array(&dim0, nullptr),
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
        accum::sparse_transform(
            SparseTransform::TO_SPARSE,
            ibmisc::make_array(&dim0, nullptr),
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

TEST_F(SpSparseTest, invert_accum)
{
    // Construct a SparseMatrix
    typedef TupleList<int, double, 2> TupleListT;
    TupleListT arr({20,10});

    arr.add({6,4}, 8.);    // Inverses of powers of 2 will be exact FP arithmetic
    arr.add({1,0}, 16.);

    TupleListT arr2;
    spcopy(accum::invert(accum::ref(arr2)), arr);

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


int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
    everytrace_init();    // Don't want everytrace for this test, it eats exceptions
    ibmisc_error = ibmisc::exception_error;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
