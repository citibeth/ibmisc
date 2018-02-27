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
#include <ibmisc/rlarray.hpp>
#include <spsparse/vector.hpp>
#include <spsparse/runlength.hpp>

using namespace ibmisc;
using namespace spsparse;
using namespace std;

static double const NaN = std::numeric_limits<double>::quiet_NaN();


// The fixture for testing class Foo.
class RunlengthTest : public ::testing::Test {
public:

    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    RunlengthTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~RunlengthTest()
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



bool eq_double(double a, double b)
{
    if (std::isnan(a)) return std::isnan(b);
    if (std::isnan(b)) return false;
    return a == b;
}


void _test_online_double(std::vector<double> const &vals)
{
    //std::vector<double> vals {1.1, 4.0, NaN, NaN, 3.0, 3.0, 3.0};

    // Online
    std::vector<int> enc_counts;
    std::vector<double> enc_values;
    {auto rle(rl_encode(
        vaccum::vector(enc_counts), vaccum::vector(enc_values), RLAlgo::PLAIN));

        for (auto v : vals) rle.add(v);
    }

    std::vector<double> vals2;
    for (auto rld(rl_decode(
        enc_counts.begin(), enc_counts.end(),
        enc_values.begin(), enc_values.end())); ++rld; )
    {
        vals2.push_back(*rld);
    }


    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        EXPECT_TRUE(eq_double(vals[i], vals2[i]));
    }

#if 0
    cout << vals << endl;
    cout << vals2 << endl;

    cout << "counts" << enc_counts << endl;
    cout << "values" << enc_values << endl;

#endif

}

void _test_online_int(std::vector<int> const &vals)
{
    // Online
    std::vector<int> enc_counts;
    std::vector<int> enc_values;
    {auto rle(rl_encode(
        vaccum::vector(enc_counts), vaccum::vector(enc_values), RLAlgo::DIFFS));

        for (auto v : vals) rle.add(v);
    }

    auto _rld(rl_decode(
        enc_counts.begin(), enc_counts.end(),
        enc_values.begin(), enc_values.end(),
        RLAlgo::DIFFS));

    std::vector<int> vals2;
    for (auto rld(rl_decode(
        enc_counts.begin(), enc_counts.end(),
        enc_values.begin(), enc_values.end(),
        RLAlgo::DIFFS)); ++rld; )
    {
        vals2.push_back(*rld);
    }


    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        EXPECT_EQ(vals[i], vals2[i]);
    }

#if 0
    cout << vals << endl;
    cout << vals2 << endl;
    cout << "counts" << enc_counts << endl;
    cout << "values" << enc_values << endl;
#endif
}


void _test_rl_vector_int(std::vector<int> const &vals)
{
    RLVector<int,int> rlv(RLAlgo::DIFFS);

    // Encode
    {auto vaccum(rlv.vaccum());
        for (auto val : vals) vaccum.add(val);
    }

    // Decode
    std::vector<int> vals2;
    for (auto gen(rlv.generator()); ++gen; ) {
        vals2.push_back(*gen);
    }

    // Compare
    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        EXPECT_EQ(vals[i], vals2[i]);
    }

#if 0
    cout << "test_rl_vector_int " << vals << endl;
    cout << "test_rl_vector_int " << vals2 << endl;
#endif
}

void _test_rl_vector_double(std::vector<double> const &vals)
{
    RLVector<int,double> rlv(RLAlgo::PLAIN);

    // Encode
    {auto vaccum(rlv.vaccum());
        for (auto val : vals) vaccum.add(val);
    }

    // Decode
    std::vector<double> vals2;
    for (auto gen(rlv.generator()); ++gen; ) {
        vals2.push_back(*gen);
    }

    // Compare
    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        EXPECT_TRUE(eq_double(vals[i], vals2[i]));
    }

#if 0
    cout << "test_rl_vector_double " << vals << endl;
    cout << "test_rl_vector_double " << vals2 << endl;
#endif
}


template<class TypeT>
void _test_rl_vector_netcdf(std::string const &fnm, RunlengthTest &self, std::vector<TypeT> const &vals, RLAlgo algo)
{
    DefaultRLEqual<TypeT> eq;

    RLVector<int,TypeT> rlv(algo);

    // Encode
    {auto vaccum(rlv.vaccum());
        for (auto val : vals) vaccum.add(val);
    }

    // Store
    std::string fname("__runlength_" + fnm + ".nc");
    self.tmpfiles.push_back(fname);
    ::remove(fname.c_str());
    {NcIO ncio(fname, 'w');
        rlv.ncio(ncio, "vals");
    }

    // Load
    RLVector<int,TypeT> rlv2;
    {NcIO ncio(fname, 'r');
        rlv2.ncio(ncio, "vals");
    }

    // Decode
    std::vector<TypeT> vals2;
    for (auto gen(rlv2.generator()); ++gen; ) {
        vals2.push_back(*gen);
    }

    // Compare
    EXPECT_EQ(vals2.size(), vals.size());
    for (int i=0; i<vals2.size(); ++i) {
        bool x = eq(vals[i], vals2[i]);
        EXPECT_TRUE(eq(vals[i], vals2[i]));
    }

#if 0
    cout << "test_rl_vector_netcdf " << vals << endl;
    cout << "test_rl_vector_netcdf " << vals2 << endl;
#endif
}


TEST_F(RunlengthTest, double)
{
    std::vector<std::vector<double>> dvalss {
        {1.1, 4.0, NaN, NaN, 3.0, 3.0, 3.0}
    };
    for (auto &vals : dvalss) {
        _test_online_double(vals);
        _test_rl_vector_double(vals);
        _test_rl_vector_netcdf<double>("double", *this, vals, RLAlgo::PLAIN);
    }
}

TEST_F(RunlengthTest, int)
{
    std::vector<std::vector<int>> ivalss {
        {1, 2, 3, 6, 7, 9, 10, 11, 17}
    };
    for (auto &vals : ivalss) {
        _test_online_int(vals);
        _test_rl_vector_int(vals);
        _test_rl_vector_netcdf<int>("int", *this, vals, RLAlgo::DIFFS);
    }
}




// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
