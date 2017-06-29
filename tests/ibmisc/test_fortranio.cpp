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
#include <iostream>
#include <cstdio>
#include <fstream>

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
        remove("sample_fortranio");
    }

//    // The mock bar library shaed by all tests
//    MockBar m_bar;
};


TEST_F(FortranIOTest, write_then_read)
{
    std::array<char, 80> str0, str1;
    blitz::Array<int, 1> vals(17);

    std::ifstream fin("sample_fortranio", ios::binary | ios::in);
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


TEST_F(FortranIOTest, cast)
{
    std::array<char, 80> str0, str1;
    blitz::Array<long, 1> vals_l(17);
    blitz::Array<int, 1> vals_i(17);

    // Read as it was written
    {std::ifstream fin("sample_fortranio", ios::binary | ios::in);
        fortran::read(fin) >> str0 >> vals_i >> fortran::endr;
    }

    // Read, casting it to long
    {std::ifstream fin("sample_fortranio", ios::binary | ios::in);
        fortran::read(fin) >> str0 >>
            fortran::blitz_cast<int, long, 1>(vals_l) >> fortran::endr;
    }

    for (int i=0; i<17; ++i) {
        EXPECT_EQ(vals_i(i), vals_l(i));
    }

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
