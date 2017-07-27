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
#include <ibmisc/bundle.hpp>

using namespace std;
using namespace ibmisc;
using namespace blitz;

// The fixture for testing class Foo.
class BundleTest : public ::testing::Test {
protected:
    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    BundleTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~BundleTest() {
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

struct ConstructAlloc {
    ArrayBundle<int,2> bundle;
    blitz::Array<int,2> a1;
    blitz::Array<int,2> a2;

    ConstructAlloc() :
        a1(bundle.add("a1", blitz::shape(2,3), {"two", "three"}, {
            std::make_tuple("longname", "Aye one"),
            std::make_tuple("units", "m")
        }, blitz::fortranArray)),
        a2(bundle.add("a2", blitz::shape(2,3), {"two", "three"}, {
            std::make_tuple("longname", "Aye two"),
            std::make_tuple("units", "m")
        }, blitz::fortranArray))
    {}

};


TEST_F(BundleTest, construct_alloc)
{
    ConstructAlloc rec;

//    EXPECT_NEQ(rec.a1.data(), rec.a2.data());
    EXPECT_EQ(1, rec.a1.lbound(0));
    EXPECT_EQ(2, rec.a1.ubound(0));
    EXPECT_EQ(1, rec.a1.lbound(1));
    EXPECT_EQ(3, rec.a1.ubound(1));
    EXPECT_EQ(6, rec.a1.size());

    std::string fname("__bundle_construct_alloc.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());

    rec.a1 = 11;
    rec.a2 = 22;

    {NcIO ncio(fname, 'w');
        rec.bundle.ncio(ncio, "", "int");
    }

    ConstructAlloc rec2;
    rec2.a1 = 17;
    {NcIO ncio(fname, 'r');
        rec2.bundle.ncio(ncio, "", "int");
    }
    EXPECT_EQ(11, rec2.a1(1,1));
    EXPECT_EQ(22, rec2.a2(1,1));
}
// -----------------------------------------------------------
struct ConstructMulti {
    ArrayBundle<int,2> bundle;
    blitz::Array<int,3> ARR;

    blitz::Array<int,2> a1;
    blitz::Array<int,2> a2;

    ConstructMulti() :
        ARR(bundle.add(blitz::shape(2,3), {"two", "three"}, {
            BundleSpec("a1", {
                make_tuple("longname", "Aye one"),
                make_tuple("units", "m")
            }),
            BundleSpec("a2", {
                make_tuple("longname", "Aye two"),
                make_tuple("units", "m")
            })
        })),
        a1(bundle.at("a1")),
        a2(bundle.at("a2"))
    {}

};

TEST_F(BundleTest, construct_multi)
{
    ConstructMulti rec;

//    EXPECT_NEQ(rec.a1.data(), rec.a2.data());
    EXPECT_EQ(0, rec.a1.lbound(0));
    EXPECT_EQ(1, rec.a1.ubound(0));
    EXPECT_EQ(0, rec.a1.lbound(1));
    EXPECT_EQ(2, rec.a1.ubound(1));
    EXPECT_EQ(6, rec.a1.size());

    std::string fname("__bundle_construct_multi.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());

    rec.a1 = 11;
    rec.a2 = 22;

    NcIO ncio(fname, 'w');
    rec.bundle.ncio(ncio, "", "int");
}

// -----------------------------------------------------------


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
