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
#include <ibmisc/iter.hpp>
#include <iostream>
#include <cstdio>
#include <memory>
#include <map>

using namespace ibmisc;

// The fixture for testing class Foo.
class IterTest : public ::testing::Test {
protected:
    // You can do set-up work for each test here.
    IterTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~IterTest() {}

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


TEST_F(IterTest, deref_iter_unique_ptr)
{
    typedef std::vector<std::unique_ptr<int>> MyVector;
    typedef DerefRandomAccessIter<int, typename MyVector::iterator> iterator;
    typedef DerefRandomAccessIter<const int, typename MyVector::const_iterator> const_iterator;

    MyVector vec1;
    vec1.push_back(std::unique_ptr<int>(new int(1)));
    vec1.push_back(std::unique_ptr<int>(new int(2)));
    vec1.push_back(std::unique_ptr<int>(new int(3)));

    int i=1;
    auto i2 = iterator(vec1.begin());
    auto i1(iterator(vec1.begin()));

    for (auto ii=iterator(vec1.begin()); ii != iterator(vec1.end()); ++ii, ++i)
        EXPECT_EQ(*ii, i);

    i=1;
    for (const_iterator ii(vec1.cbegin()); ii != const_iterator(vec1.cend()); ++ii, ++i)
        EXPECT_EQ(*ii, i);

}
// -----------------------------------------------------------
TEST_F(IterTest, deref_iter_cptr)
{
    std::vector<int> ints(3);
    for (int i=0; i<3; ++i) ints[i] = i+1;

    typedef std::vector<int *> MyVector;
    typedef DerefRandomAccessIter<int, typename MyVector::iterator> iterator;
    typedef DerefRandomAccessIter<const int, typename MyVector::const_iterator> const_iterator;

    MyVector vec1;
    vec1.push_back(&ints[0]);
    vec1.push_back(&ints[1]);
    vec1.push_back(&ints[2]);

    int i=1;
    for (auto ii=iterator(vec1.begin()); ii != iterator(vec1.end()); ++ii, ++i)
        EXPECT_EQ(*ii, i);
}
// -----------------------------------------------------------
TEST_F(IterTest, second_iter)
{
    typedef std::map<int, double> MyMap;

    MyMap map1;
    map1.insert(std::make_pair(1, 1.0));
    map1.insert(std::make_pair(2, 2.0));
    map1.insert(std::make_pair(3, 3.0));

    typedef SecondIter<int, double, typename MyMap::iterator> iterator;
    typedef SecondIter<const int, const double, typename MyMap::const_iterator> const_iterator;

    int i=1;
    for (auto ii=iterator(map1.begin()); ii != iterator(map1.end()); ++ii, ++i) {
        EXPECT_EQ(*ii, (double)i);
        EXPECT_EQ(ii.key(), i);
    }
}
// -----------------------------------------------------------


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
