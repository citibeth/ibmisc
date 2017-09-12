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

#include <ibmisc/Test.hpp>
#include <ibmisc/permutation.hpp>
#include <iostream>
#include <memory>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <functional>

using namespace ibmisc;

// The fixture for testing class Foo.
class PermutationTest : public ibmisc::Test {
protected:

    // You can do set-up work for each test here.
    PermutationTest() : ibmisc::Test(false) {}    // keep=true
};

TEST_F(PermutationTest, sort_permutation)
{

    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine(rnd_device());
    std::uniform_int_distribution<int> dist(1, 52);

    auto gen = std::bind(dist, mersenne_engine);

    size_t N=100;
    std::vector<int> vec(N);
    std::generate(std::begin(vec), std::end(vec), gen);


    // Optional

    std::vector<size_t> sorted_perm;
    ibmisc::sorted_permutation(std::begin(vec), std::end(vec), sorted_perm);


    auto ii(vec.begin());
    int last = vec[sorted_perm[0]];
    for (int i=1; i<N; ++i) {
        EXPECT_GE(vec[sorted_perm[i]], last);
        last = vec[sorted_perm[i]];
    }
}

// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
