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
#include <ibmisc/fortranio.hpp>
#include <cstdio>

using namespace std;
using namespace ibmisc;

extern "C" void write_sample();    // Defined in help_fortranio.F90

// The fixture for testing class Foo.
class FortranIOTest : public ::testing::Test {
protected:
    // You can do set-up work for each test here.
    FortranIOTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~FortranIOTest() {}

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() {
        write_sample();
    }

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown() {
        remove("sample_fortranio_le");
        remove("sample_fortranio_be");
    }

//    // The mock bar library shaed by all tests
//    MockBar m_bar;
};




void test_write_then_read(std::string const &fname, Endian endian)
{
    std::array<char, 80> str0, str1;
    blitz::Array<int, 1> vals(17);

    fortran::UnformattedInput fin(fname, endian);
    fortran::read(fin) >> str0 >> vals >> fortran::endr;

    EXPECT_EQ("Hello World", fortran::trim(str0));
    for (int i=0; i<17; ++i) {
        EXPECT_EQ(i+1, vals(i));
    }

    fortran::read(fin) >> fortran::endr;    // skip

    fortran::read(fin) >> str1 >> vals >> str0 >> fortran::endr;
    EXPECT_EQ("Hello World", fortran::trim(str1));
    for (int i=0; i<17; ++i) {
        EXPECT_EQ(17-(i+1), vals(i));
    }
    EXPECT_EQ("Goodbye World", fortran::trim(str0));
    
}

void test_cast(std::string const &fname, ibmisc::Endian endian)
{
    std::array<char, 80> str0, str1;
    blitz::Array<long, 1> vals_l(17);
    blitz::Array<int, 1> vals_i(17);

    // Read as it was written
    {fortran::UnformattedInput fin(fname, endian);
        fortran::read(fin) >> str0 >> vals_i >> fortran::endr;
    }

    // Read, casting it to long
    {fortran::UnformattedInput fin(fname, endian);
        fortran::read(fin) >> str0 >>
            fortran::blitz_cast<int, long, 1>(vals_l) >> fortran::endr;
    }

    for (int i=0; i<17; ++i) {
        EXPECT_EQ(vals_i(i), vals_l(i));
    }

}

void test_shape_cast(std::string const &fname, ibmisc::Endian endian)
{
    using namespace fortran;

    std::array<char, 80> str0, str1;
    std::vector<Shape<2>> stdshapes {
        Shape<2>({"im1", "jm1"}, {1,4}),
        Shape<2>({"im2", "jm2"}, {4,3})
    };

    // Don't retrieve shape
    {fortran::UnformattedInput fin(fname, endian);
        // Skip 3 records
        for (int i=0; i<3; ++i)
            fortran::read(fin) >> fortran::endr;

        // Read record [4] using wildcard
        blitz::Array<int,2> data;
        fortran::read(fin)
            >> str0
            >> fortran::allocate<int,2>(data, stdshapes)
            >> str1
            >> fortran::endr;
        EXPECT_EQ(1, data.lbound(0));
        EXPECT_EQ(4, data.ubound(0));
        EXPECT_EQ(1, data.lbound(1));
        EXPECT_EQ(3, data.ubound(1));
    }



    // Retrieve shape
    {fortran::UnformattedInput fin(fname, endian);
        // Skip 3 records
        for (int i=0; i<3; ++i)
            fortran::read(fin) >> fortran::endr;

        // Read record [4] using wildcard
        blitz::Array<int,2> data;
        Shape<2> const *data_shape;   // Tell us what it did
        fortran::read(fin)
            >> str0
            >> fortran::allocate<int,2>(data, data_shape, stdshapes)
            >> str1
            >> fortran::endr;

        EXPECT_EQ(1, data.lbound(0));
        EXPECT_EQ(4, data.ubound(0));
        EXPECT_EQ(1, data.lbound(1));
        EXPECT_EQ(3, data.ubound(1));
        for (int i=1; i<=3; ++i) {
        for (int j=1; j<=3; ++j) {
            EXPECT_EQ((i-1)*j, data(i,j));
        }}
        EXPECT_EQ("im2", data_shape->sshape[0]);
        EXPECT_EQ("jm2", data_shape->sshape[1]);
    }





}


TEST_F(FortranIOTest, shape_cast_be)
    { test_shape_cast("sample_fortranio_be", ibmisc::Endian::BIG); }
TEST_F(FortranIOTest, shape_cast_le)
    { test_shape_cast("sample_fortranio_le", ibmisc::Endian::LITTLE); }







TEST_F(FortranIOTest, write_then_read_be)
{
    test_write_then_read("sample_fortranio_be", ibmisc::Endian::BIG);

}
TEST_F(FortranIOTest, write_then_read_le)
{
    test_write_then_read("sample_fortranio_le", ibmisc::Endian::LITTLE);
}


TEST_F(FortranIOTest, cast_be)
{
    test_cast("sample_fortranio_be", ibmisc::Endian::BIG);
}
TEST_F(FortranIOTest, cast_le)
{
    test_cast("sample_fortranio_le", ibmisc::Endian::LITTLE);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
