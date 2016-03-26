/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/blitz.hpp>
#include <iostream>
#include <cstdio>

using namespace ibmisc;
using namespace blitz;

// The fixture for testing class Foo.
class BlitzTest : public ::testing::Test {
protected:

    // You can do set-up work for each test here.
    BlitzTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~BlitzTest()
    {
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    virtual void SetUp() {}

    // Code here will be called immediately after each test (right
    // before the destructor).
    virtual void TearDown() {}

};

TEST_F(BlitzTest, tiny_conv)
{
    std::array<int, 2> aval = {0, 1};
    auto bval(ibmisc::to_tiny<long, int,2>(aval));

    blitz::TinyVector<long,2> bval2;
    ibmisc::to_tiny<long, int, 2>(bval2, aval);
}

TEST_F(BlitzTest, blitz_conv)
{
    std::vector<int> v = {1, 2, 3};
    auto b(to_blitz(v));
    EXPECT_EQ(b.extent(0), v.size());
    EXPECT_EQ(&b(0), &v[0]);

    std::vector<double> v2 = {1, 2, 3};
    auto b2(to_blitz(v));
    EXPECT_EQ(b.extent(0), v.size());
    EXPECT_EQ(&b(0), &v[0]);
    EXPECT_EQ(&b(1), &v[1]);
    EXPECT_EQ(&b(2), &v[2]);

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
