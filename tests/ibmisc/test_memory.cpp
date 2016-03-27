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
#include <ibmisc/memory.hpp>
#include <ibmisc/string.hpp>
#include <iostream>
#include <cstdio>
#include <memory>
#include <map>

using namespace std;
using namespace ibmisc;

// The fixture for testing class Foo.
class MemoryTest : public ::testing::Test {
protected:
    // You can do set-up work for each test here.
    MemoryTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~MemoryTest() {}

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

string hello = "hello";
string *get_hello()
{
    return &hello; 
}
unique_ptr<string> new_hello()
    { return unique_ptr<string>(new string("hello")); }

TEST_F(MemoryTest, lazy_ptr_test)
{
    LazyPtr<string> p0(get_hello());
    LazyPtr<string> p1(new_hello());
    LazyPtr<string> p2(&get_hello);
    LazyPtr<string> p3(&new_hello);

    EXPECT_EQ("hello", *p0);
    EXPECT_EQ("hello", *p1);
    EXPECT_EQ("hello", *p2);
    EXPECT_EQ("hello", *p3);
}

TEST_F(MemoryTest, string_printf)
{
    std::string s1(string_printf("Hello World"));
    EXPECT_EQ("Hello World", s1);

    std::string s2(string_printf("H%s%d", "ello ", 17));
    EXPECT_EQ("Hello 17", s2);
}


// -----------------------------------------------------------
// -----------------------------------------------------------


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
