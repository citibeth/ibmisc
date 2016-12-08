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
#include <ibmisc/indexing.hpp>
#include <ibmisc/IndexSet.hpp>
#include <ibmisc/Test.hpp>
#include <iostream>
#include <memory>

using namespace ibmisc;
using namespace netCDF;

// The fixture for testing class Foo.
class IndexingTest : public ibmisc::Test {
protected:

    // You can do set-up work for each test here.
    IndexingTest() : ibmisc::Test(true) {}    // keep=true
};

TEST_F(IndexingTest, indexing_column_major_test)
{
    Indexing ind(
        {"d0", "d1"},    // Names
        {0,0},      // Base
        {5,4},      // Extent
        {1,0});     // Column major
    std::array<int,2> tuple{3,2};
    long ix2 = ind.tuple_to_index(tuple);
    auto tuple2(ind.index_to_tuple<int,2>(ix2));

    EXPECT_EQ(20, ind.extent());
    EXPECT_EQ(13, ix2);
    EXPECT_EQ(tuple, tuple2);

}

TEST_F(IndexingTest, indexing_row_major_test)
{
    Indexing ind(
        {"d0", "d1"},
        {0,0},      // Base
        {4,5},      // Extent
        {0,1});     // Row major
    std::array<int,2> tuple;

    tuple = {3,2};
    long ix2 = ind.tuple_to_index(tuple);
    auto tuple2(ind.index_to_tuple<int,2>(ix2));

    EXPECT_EQ(20, ind.extent());
    EXPECT_EQ(17, ix2);
    EXPECT_EQ(tuple, tuple2);

}

TEST_F(IndexingTest, indexing_netcdf)
{
    std::string fname(tmp_fname("__netcdf_indexing_test.nc"));


    Indexing ind(
        {"d0", "d1"},
        {0,0},      // Base
        {4,5},      // Extent
        {0,1});     // Row major
    {
        NcIO ncio(fname, NcFile::replace);
        ind.ncio(ncio, "indexing");  
        ncio.close();
    }

    Indexing ind2;
    {
        NcIO ncio(fname, 'r');
        ind2.ncio(ncio, "indexing");
        ncio.close();
    }

    EXPECT_TRUE(ind2 == ind);

}

TEST_F(IndexingTest, make_blitzcol_major)
{
    Indexing ind(
        {"d0", "d1"},
        {1,2},      // Base
        {4,5},      // Extent
        {1,0});     // Column major

    auto arr(ind.make_blitz<double,2>());
    EXPECT_EQ(arr.lbound(0), 1);
    EXPECT_EQ(arr.lbound(1), 2);
    EXPECT_EQ(arr.ubound(0), 4);
    EXPECT_EQ(arr.ubound(1), 6);
    EXPECT_EQ(1, &arr(2,2) - &arr(1,2));
    EXPECT_EQ(4, &arr(1,3) - &arr(1,2));

}
TEST_F(IndexingTest, make_blitzrow_major)
{
    Indexing ind(
        {"d0", "d1"},
        {1,2},      // Base
        {4,5},      // Extent
        {0,1});     // Column major

    auto arr(ind.make_blitz<double,2>());
    EXPECT_EQ(arr.lbound(0), 1);
    EXPECT_EQ(arr.lbound(1), 2);
    EXPECT_EQ(arr.ubound(0), 4);
    EXPECT_EQ(arr.ubound(1), 6);
    EXPECT_EQ(5, &arr(2,2) - &arr(1,2));
    EXPECT_EQ(1, &arr(1,3) - &arr(1,2));
}

TEST_F(IndexingTest, to_blitzcol_major)
{
    Indexing ind(
        {"d0", "d1"},
        {1,2},      // Base
        {4,5},      // Extent
        {1,0});     // Column major

    double mem[20];
    auto arr(ind.to_blitz<double,2>(&mem[0]));
    EXPECT_EQ(arr.lbound(0), 1);
    EXPECT_EQ(arr.lbound(1), 2);
    EXPECT_EQ(arr.ubound(0), 4);
    EXPECT_EQ(arr.ubound(1), 6);
    EXPECT_EQ(1, &arr(2,2) - &arr(1,2));
    EXPECT_EQ(4, &arr(1,3) - &arr(1,2));

}
TEST_F(IndexingTest, to_blitzrow_major)
{
    Indexing ind(
        {"d0", "d1"},
        {1,2},      // Base
        {4,5},      // Extent
        {0,1});     // Column major

    double mem[20];
    auto arr(ind.to_blitz<double,2>(&mem[0]));
    EXPECT_EQ(arr.lbound(0), 1);
    EXPECT_EQ(arr.lbound(1), 2);
    EXPECT_EQ(arr.ubound(0), 4);
    EXPECT_EQ(arr.ubound(1), 6);
    EXPECT_EQ(5, &arr(2,2) - &arr(1,2));
    EXPECT_EQ(1, &arr(1,3) - &arr(1,2));

}

// -----------------------------------------------------------
TEST_F(IndexingTest, domain)
{
    Indexing ind(
        {"d0", "d1"},
        {0,0},      // Base
        {4,5},      // Extent
        {0,1});     // Row major

    Domain domain({3,3}, {4,5});

    bool x;
    x = domain.in_domain<int,2>({3,4}); EXPECT_TRUE(x);
    x = domain.in_domain<int,2>({3,5}); EXPECT_FALSE(x);
    x = domain.in_domain<int,2>({0,0}); EXPECT_FALSE(x);

    x = in_domain(&domain, &ind, 0L); EXPECT_FALSE(x);
    x = in_domain(&domain, &ind, 15L); EXPECT_FALSE(x);
    x = in_domain(&domain, &ind, 19L); EXPECT_TRUE(x);
    x = in_domain(&domain, &ind, 20L); EXPECT_FALSE(x);

}

TEST_F(IndexingTest, domain_netcdf)
{
    std::string fname(tmp_fname("__netcdf_domain_test.nc"));

    Domain domain({3,3}, {4,5});

    {
        NcIO ncio(fname, NcFile::replace);
        domain.ncio(ncio, "domain");  
        ncio.close();
    }

    Domain domain2;
    {
        NcIO ncio(fname, 'r');
        domain2.ncio(ncio, "domain");
        ncio.close();
    }

    EXPECT_TRUE(domain2 == domain);

}


// -----------------------------------------------------------


TEST_F(IndexingTest, index_set)
{
    IndexSet<std::string> names;

    names.insert("A");
    names.insert("B");
    names.insert("C");

    EXPECT_EQ("A", names[0]);
    EXPECT_EQ("B", names[1]);
    EXPECT_EQ("C", names[2]);

    EXPECT_EQ(0, names.at("A"));
    EXPECT_EQ(1, names.at("B"));
    EXPECT_EQ(2, names.at("C"));

    auto ii(names.begin());
    EXPECT_EQ(names[0], *ii);
    ++ii;
    EXPECT_EQ(names[1], *ii);
    ++ii;
    EXPECT_EQ(names[2], *ii);
    ++ii;
    EXPECT_EQ(names.end(), ii);

    EXPECT_TRUE(names.contains("A"));
    EXPECT_FALSE(names.contains("ZZTOP"));
    EXPECT_THROW(names.at("ZZTOP"), ibmisc::Exception);
    EXPECT_THROW(names[17], ibmisc::Exception);
}


// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
