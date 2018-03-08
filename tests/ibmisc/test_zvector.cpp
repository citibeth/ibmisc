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

#include <gtest/gtest.h>
#include <ibmisc/Test.hpp>
#include <ibmisc/array.hpp>
#include <prettyprint.hpp>
#include <iostream>
#include <ibmisc/blitz.hpp>
//#include <ibmisc/rlarray.hpp>
#include <spsparse/zvector.hpp>
#include <spsparse/eigen.hpp>
#include <ibmisc/zarray.hpp>



// -----------------------------------------------------------------
/** Equality tester takes NaN into account for double and float.
Use std::equal_to<TypeT> instead, if you don't like this. */
template<class TypeT>
struct DefaultZVEqual {
    bool operator()(TypeT const &a, TypeT const &b) const
        { return (a == b); }
};

template<>
bool DefaultZVEqual<double>::operator()(double const &a, double const &b) const
{
    switch(
        (std::isnan(a) ? 2 : 0) +
        (std::isnan(b) ? 1 : 0))
    {
        case 3:
            return true;
        case 1:
        case 2:
            return false;
        case 0:
            return (a == b);
    }
}

template<>
bool DefaultZVEqual<float>::operator()(float const &a, float const &b) const
{
    switch(
        (std::isnan(a) ? 2 : 0) +
        (std::isnan(b) ? 1 : 0))
    {
        case 3:
            return true;
        case 1:
        case 2:
            return false;
        case 0:
            return (a == b);
    }
}
// -----------------------------------------------------------------



using namespace ibmisc;
using namespace spsparse;
using namespace std;

static double const NaN = std::numeric_limits<double>::quiet_NaN();


// The fixture for testing class Foo.
class ZVectorTest : public ::testing::Test {
public:

    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    ZVectorTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~ZVectorTest()
    {
        for (auto ii(tmpfiles.begin()); ii != tmpfiles.end(); ++ii) {
//          ::remove(ii->c_str());
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


template<class TypeT,int RANK>
void _test_zvector(std::vector<std::array<TypeT,RANK>> const &vals, ZVAlgo algo)
{
    std::vector<char> zbuf;

    // Encode
    {vaccum::ZVector<TypeT,RANK> accum(zbuf, algo);
        for (auto val : vals) accum.add(val);
    }

    printf("Encoded to %ld bytes\n", zbuf.size());

    // Decode
    std::vector<std::array<TypeT,RANK>> vals2;
    for (vgen::ZVector<TypeT,RANK> gen(zbuf); ++gen; ) {
        vals2.push_back(*gen);
    }

    // Compare
    DefaultZVEqual<TypeT> eq;
    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        for (int j=0; j<RANK; ++j) {
            EXPECT_TRUE(eq(vals[i][j], vals2[i][j]));
        }
    }

#if 1
    cout << "test_zvector vals  " << vals << endl;
    cout << "test_zvector vals2 " << vals2 << endl;
#endif
}



TEST_F(ZVectorTest, int)
{
    std::vector<std::vector<std::array<int,2>>> ivalss {
        {},
        {{0,1}, {0,2}, {0,3}, {0,6}, {1,7}, {1,9}, {1,10}, {1,11}, {1,17}}
    };
    for (auto &vals : ivalss) {
        _test_zvector<int,2>(vals, ZVAlgo::PLAIN);
        _test_zvector<int,2>(vals, ZVAlgo::DIFFS);
    }
}

TEST_F(ZVectorTest, double)
{
    std::vector<std::vector<std::array<double,1>>> dvalss {
        {},
        {{1.1}, {4.0}, {NaN}, {NaN}, {3.0}, {3.0}, {3.0}}
    };
    for (auto &vals : dvalss) {
        _test_zvector<double,1>(vals, ZVAlgo::PLAIN);
    }
}
// -----------------------------------------------------------------

#if 1
TEST_F(ZVectorTest, ZArray1)
{
    // Set up a TupleList
    TupleList<int, double, 2> arr1({400,400});
    arr1.add({0,1}, 1.);
    arr1.add({0,2}, 1.);
    arr1.add({0,3}, 1.);
    arr1.add({0,4}, .5);
    arr1.add({1,9}, .6);
    arr1.add({1,10}, 1.1);
    arr1.add({1,11}, 1.1);
    arr1.add({1,12}, 1.1);
    arr1.add({1,13}, 1.0);
    arr1.add({2,25}, 1.9);
    arr1.add({2,26}, 1.9);
    arr1.add({2,27}, 1.9);

    // ZVector-Compress it
    ZArray<int,double,2> zsa1(arr1.shape());
    {auto accum(zsa1.accum());
        for (auto &tup : arr1) accum.add(tup.index(), tup.value());
    }

    // Store it
    std::string fname("__zsparsearray1.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());
    {NcIO ncio(fname, 'w');
        zsa1.ncio(ncio, "vals");
    }

    // Load it
    ZArray<int,double,2> zsa2;
    {NcIO ncio(fname, 'r');
        zsa2.ncio(ncio, "vals");
    }
    EXPECT_EQ(zsa1.nnz(), zsa2.nnz());
    EXPECT_EQ(zsa1.shape()[0], zsa2.shape()[0]);
    EXPECT_EQ(zsa1.shape()[1], zsa2.shape()[1]);


    // Uncompress it
    TupleList<int,double,2> arr2(zsa2.shape());
    for (auto ii(zsa2.generator()); ++*ii; ) {
        arr2.add(ii->index(), ii->value());
    }

    // Compare
    EXPECT_EQ(arr1.shape(), arr2.shape());
    for (size_t i=0; i<arr1.size(); ++i) {
        Tuple<int,double,2> const &tup1(arr1.tuples[i]);
        Tuple<int,double,2> const &tup2(arr2.tuples[i]);
        bool const eq = (tup1 == tup2);
        EXPECT_TRUE(eq);
    }
}
#endif


// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
