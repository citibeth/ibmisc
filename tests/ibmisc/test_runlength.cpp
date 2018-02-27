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
#include <ibmisc/runlength.hpp>
#include <ibmisc/array.hpp>
#include <prettyprint.hpp>
#include <iostream>
#include <ibmisc/blitz.hpp>
#include <spsparse/vector.hpp>

using namespace ibmisc;
using namespace spsparse;
using namespace std;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

// The fixture for testing class Foo.
class RunlengthTest : public ibmisc::Test {
protected:

    // You can do set-up work for each test here.
    RunlengthTest() : ibmisc::Test(false) {}    // keep=true
};

bool eq_double(double a, double b)
{
    if (std::isnan(a)) return std::isnan(b);
    if (std::isnan(b)) return false;
    return a == b;
}


void _test_sample_double(std::vector<double> const &vals)
{
    // std::vector<double> vals {1.1, 4.0, NaN, NaN, 3.0, 3.0, 3.0};
    EqualUsingNaN eq;
    auto enc(rlencode(to_blitz(vals), false, eq));
    auto vals2(rldecode(to_blitz(enc.ends), to_blitz(enc.values)));

    EXPECT_EQ(vals2.extent(0), vals.size());
    for (int i=0; i<vals2.extent(0); ++i) {
        EXPECT_TRUE(eq_double(vals[i], vals2(i)));
    }
}

void _test_sample_int(std::vector<int> const &vals)
{
    // std::vector<int> vals {1, 2, 3, 6, 7, 9, 10, 11, 17};
    auto enc(rlencode(to_blitz(vals), true));

#if 0
    cout << enc.ends << endl;
    cout << enc.values << endl;
#endif

    auto vals2(rldecode(to_blitz(enc.ends), to_blitz(enc.values), true));

    EXPECT_EQ(vals2.extent(0), vals.size());
    for (int i=0; i<vals2.extent(0); ++i) {
        EXPECT_EQ(vals[i], vals2(i));
    }
}

void _test_online_double(std::vector<double> const &vals)
{
    //std::vector<double> vals {1.1, 4.0, NaN, NaN, 3.0, 3.0, 3.0};
    EqualUsingNaN eq;

    // Online
    std::vector<int> enc_counts;
    std::vector<double> enc_values;
    {auto rle(rl_encoder(
        accum::vector(enc_counts), accum::vector(enc_values), false, EqualUsingNaN()));

        for (auto v : vals) rle.add(v);
    }

    auto _rld(rl_decoder(
        enc_counts.begin(), enc_counts.end(),
        enc_values.begin(), enc_values.end()));

    std::vector<double> vals2;
    for (auto rld(std::move(_rld)); ++rld; ) {
        vals2.push_back(*rld);
    }



    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        EXPECT_TRUE(eq_double(vals[i], vals2[i]));
    }

#if 0
    cout << vals << endl;
    cout << vals2 << endl;
#endif

}

void _test_online_int(std::vector<int> const &vals)
{
    // Online
    std::vector<int> enc_counts;
    std::vector<int> enc_values;
    {auto rle(rl_encoder(
        accum::vector(enc_counts), accum::vector(enc_values), true));

        for (auto v : vals) rle.add(v);
    }

    auto _rld(rl_decoder(
        enc_counts.begin(), enc_counts.end(),
        enc_values.begin(), enc_values.end(),
        true));

    std::vector<int> vals2;
    for (auto rld(std::move(_rld)); ++rld; ) {
        vals2.push_back(*rld);
    }



    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        EXPECT_EQ(vals[i], vals2[i]);
    }

#if 0
    cout << vals << endl;
    cout << vals2 << endl;
#endif
}


TEST_F(RunlengthTest, double)
{
    std::vector<std::vector<double>> dvalss {
        {1.1, 4.0, NaN, NaN, 3.0, 3.0, 3.0}
    };
    for (auto &vals : dvalss) {
        _test_sample_double(vals);
        _test_online_double(vals);
    }
}

TEST_F(RunlengthTest, int)
{
    std::vector<std::vector<int>> ivalss {
        {1, 2, 3, 6, 7, 9, 10, 11, 17}
    };
    for (auto &vals : ivalss) {
        _test_sample_int(vals);
        _test_online_int(vals);
    }
}




// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
