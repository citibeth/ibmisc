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
#include <ibmisc/netcdf.hpp>
#include <iostream>
#include <cstdio>
#include <netcdf>
#include <everytrace.h>

using namespace ibmisc;
using namespace netCDF;

// The fixture for testing class Foo.
class NetcdfTest : public ::testing::Test {
protected:

    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    NetcdfTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~NetcdfTest()
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

TEST_F(NetcdfTest, get_or_add_dim) {

    // If the constructor and destructor are not enough for setting up
    std::string fname("__netcdf_get_or_add_dim_test.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());
    NcIO ncio(fname, 'w');
    NcDim dim1 = get_or_add_dim(ncio, "dim1", 1);
    NcDim dim2 = get_or_add_dim(ncio, "dim2", 2);
    NcDim dimU = get_or_add_dim(ncio, "dimU");

    EXPECT_EQ(1, dim1.getSize());
    EXPECT_EQ(2, dim2.getSize());
    EXPECT_FALSE(dim1.isUnlimited());
    EXPECT_TRUE(dimU.isUnlimited());

    NcDim dim1b = get_or_add_dim(ncio, "dim1", 1);
    EXPECT_EQ(1, dim1b.getSize());
    EXPECT_FALSE(dim1b.isUnlimited());

    NcDim dimUb = get_or_add_dim(ncio, "dimU");
    EXPECT_TRUE(dimUb.isUnlimited());

    // Try to reset to unlimited
    EXPECT_THROW(get_or_add_dim(ncio, "dim1"), ibmisc::Exception);
    // Try to reset to another value
    EXPECT_THROW(get_or_add_dim(ncio, "dim1", 17), ibmisc::Exception);
    // Try to reset from unlimited
    EXPECT_THROW(get_or_add_dim(ncio, "dimU", 17), ibmisc::Exception);

    get_or_add_var(ncio, "var1", "int", {dim1});
    get_or_add_var(ncio, "var1", "int", {dim1});
    get_or_add_var(ncio, "var1", "int", {dim1b});
    NcDim dim1c = get_or_add_dim(ncio, "dim1c", 1);
    EXPECT_THROW(get_or_add_var(ncio, "var1", "int", {dim1c}), ibmisc::Exception);
    EXPECT_THROW(get_or_add_var(ncio, "var1", "int", {dim2}), ibmisc::Exception);

}

TEST_F(NetcdfTest, blitz)
{
    std::string fname("__netcdf_blitz_test.nc");
    tmpfiles.push_back(fname);

    ::remove(fname.c_str());
    
    blitz::Array<double,2> A(4,5);
    for (int i=0; i<A.extent(0); ++i) {
    for (int j=0; j<A.extent(1); ++j) {
      A(i,j) = i*j;
    }}

    blitz::Array<double,2> B(A*2);
    blitz::Array<double,2> B_f(c_to_f(B));

    blitz::Array<double,1> C(17);
    for (int i=0; i<C.extent(0); ++i) C(i) = 52 + i*2;
    blitz::Array<double,1> C_f(c_to_f(C));

    std::vector<std::string> strings = {"s1", "s2"};

    // ---------- Write
    {
    ibmisc::NcIO ncio(fname, 'w');
    auto dims = ibmisc::get_or_add_dims(ncio, A, {"dim4", "dim5"});
    ibmisc::ncio_blitz(ncio, A, true, "A", "double", dims);
    auto dims_f = ibmisc::get_or_add_dims(ncio, B_f, {"dim5", "dim4"});
    ibmisc::ncio_blitz(ncio, B_f, true, "B", "double", dims_f);

    auto dimsC = ibmisc::get_or_add_dims(ncio, C, {"dim17"});
    ibmisc::ncio_blitz(ncio, C, true, "C", "double", dimsC);

    auto info_v = get_or_add_var(ncio, "info", "int64", {});
    get_or_put_att(info_v, ncio.rw, "strings", "", strings);    // Type ignored here

    ncio.close();
    }

    // ---------- Read
    ibmisc::NcIO ncio(fname, NcFile::read);

    blitz::Array<double,2> A2;
    blitz::Array<double,1> C2;
    std::vector<std::string> strings2;
    auto dims = ibmisc::get_or_add_dims(ncio, A, {"dim4", "dim5"});
    ibmisc::ncio_blitz(ncio, A2, true, "A", "double", dims);
    auto B2(ibmisc::nc_read_blitz<double,2>(ncio.nc, "B"));
//    auto dimsC = ibmisc::get_or_add_dims(ncio, C, {"dim17"});
    auto dimsC = ibmisc::get_dims(ncio, {"dim17"});
    ibmisc::ncio_blitz(ncio, C2, true, "C", "double", dimsC);

    auto info_v = get_or_add_var(ncio, "info", "int64", {});
    get_or_put_att(info_v, ncio.rw, "strings", "", strings2);

    ncio.close();
    EXPECT_EQ(A.extent(0), A2.extent(0));
    for (int i=0; i<A.extent(0); ++i) {
    for (int j=0; j<A.extent(1); ++j) {
        EXPECT_EQ(A(i,j), A2(i,j));
    }}
    EXPECT_EQ(strings.size(), strings2.size());
    for (size_t i=0; i<strings.size(); ++i) {
        EXPECT_EQ(strings[i], strings2[i]);
    }

    EXPECT_EQ(C.extent(0), C2.extent(0));
    for (int i=0; i<C.extent(0); ++i) {
        EXPECT_EQ(C(i), C2(i));
    }

}

TEST_F(NetcdfTest, vector)
{
    std::string fname("__netcdf_vector_test.nc");
    tmpfiles.push_back(fname);

    ::remove(fname.c_str());

    std::vector<double> vec;
    for (int i=0; i<4; ++i) vec.push_back(i+1);

    // ---------- Write
    {
        ibmisc::NcIO ncio(fname, 'w');
        auto dim = ibmisc::get_or_add_dim(ncio, "dim1", vec.size());
        ibmisc::ncio_vector(ncio, vec, true, "vec", "double", {dim});
        ncio.close();
    }

    {
        std::vector<double> vec2;
        ibmisc::NcIO ncio(fname, NcFile::read);
        auto dim = ibmisc::get_or_add_dim(ncio, "dim1", vec2.size());
        ibmisc::ncio_vector(ncio, vec2, true, "vec", "double", {dim});
        ncio.close();

        EXPECT_EQ(vec, vec2);
    }

}


int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
    everytrace_init();    // Don't want everytrace for this test, it eats exceptions
    ibmisc_error = ibmisc::exception_error;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
