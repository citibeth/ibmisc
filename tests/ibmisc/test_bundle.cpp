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



struct MyClass_Bundle : public ArrayBundle<int,2> {
    MyClass_Bundle() : ArrayBundle<int,2>({
        def("a1", {2,3}, {"two", "three"}, {
            "longname", "Aye one",
            "units", "m"
        }),
        def("a2", {2,3}, {"two", "three"}, {
            "longname", "Aye two",
            "units", "m"
        }),

    }) {}
};

struct MyClass {
    MyClass_Bundle bundle;
    blitz::Array<int,2> &a1;
    blitz::Array<int,2> &a2;

    MyClass() :
        a1(bundle.array("a1")),
        a2(bundle.array("a2"))
    {
        // bundle.allocate(true, blitz::fortranArray);
    }
};


TEST_F(BundleTest, construct_alloc)
{
    // Allocate as C (row-major)
    MyClass rec_c;
    rec_c.bundle.allocate(true);
    EXPECT_EQ(0, rec_c.a1.lbound(0));
    EXPECT_EQ(1, rec_c.a1.ubound(0));
    EXPECT_EQ(0, rec_c.a1.lbound(1));
    EXPECT_EQ(2, rec_c.a1.ubound(1));
    EXPECT_EQ(6, rec_c.a1.size());
    EXPECT_TRUE(rec_c.a1.stride(0) > rec_c.a1.stride(1));

    // Allocate as Fortran (column-major)
    MyClass rec_f;
    rec_f.bundle.allocate(true, fortranArray);
    EXPECT_EQ(1, rec_f.a1.lbound(0));
    EXPECT_EQ(2, rec_f.a1.ubound(0));
    EXPECT_EQ(1, rec_f.a1.lbound(1));
    EXPECT_EQ(3, rec_f.a1.ubound(1));
    EXPECT_EQ(6, rec_f.a1.size());
    EXPECT_TRUE(rec_f.a1.stride(1) > rec_f.a1.stride(0));



    std::string fname("__bundle_construct_alloc.nc");
//    tmpfiles.push_back(fname);
    ::remove(fname.c_str());

    rec_f.a1 = 11;
    rec_f.a2 = 22;

    // Write it out
    {NcIO ncio(fname, 'w');
//        rec_f.bundle.ncio(ncio, {}, false, "", "int");
//        using namespace std::placeholders;
//        rec_f.bundle.ncio(ncio, {}, "", "int",
//            std::bind(&_ncio_blitz::_whole1, _1, _2, _3,
//            std::vector<netCDF::NcDim>{},false,true));
        rec_f.bundle.ncio_whole(ncio, {}, "", "int");
    }

#if 0
    // Read it back, allocating as we go (as fortranArray)
    {
        MyClass rec2;
        rec2.a1 = 17;
        {NcIO ncio(fname, 'r');
            rec2.bundle.ncio(ncio, {}, true, "", "int", fortranArray);
        }
        EXPECT_EQ(1, rec2.a1.lbound(0));
        EXPECT_EQ(2, rec2.a1.ubound(0));
        EXPECT_EQ(1, rec2.a1.lbound(1));
        EXPECT_EQ(3, rec2.a1.ubound(1));
        EXPECT_EQ(6, rec2.a1.size());
        EXPECT_EQ(11, rec2.a1(1,1));
        EXPECT_EQ(22, rec2.a2(1,1));
    }

    // Read it back, allocating as we go (as C array; indices are reversed)
    {
        MyClass rec2;
        rec2.a1 = 17;
        {NcIO ncio(fname, 'r');
            rec2.bundle.ncio(ncio, {}, true, "", "int");
        }
        EXPECT_EQ(0, rec2.a1.lbound(1));
        EXPECT_EQ(1, rec2.a1.ubound(1));
        EXPECT_EQ(0, rec2.a1.lbound(0));
        EXPECT_EQ(2, rec2.a1.ubound(0));
        EXPECT_EQ(6, rec2.a1.size());
        EXPECT_EQ(11, rec2.a1(1,1));
        EXPECT_EQ(22, rec2.a2(1,1));
    }

#endif
}
#if 0
// -----------------------------------------------------------
struct MyClass2_Bundle : public ArrayBundle<int,2> {
    MyClass2_Bundle() : ArrayBundle<int,2>({
        def("a1", {
            "longname", "Aye one",
            "units", "m"
        }),
        def("a2", {
            "longname", "Aye two",
            "units", "m"
        }),

    }) {}
};


TEST_F(BundleTest, late_shape)
{
    MyClass2_Bundle bundle;
    bundle.at("a1").allocate(blitz::shape(2,3), {"two", "three"});
    bundle.at("a2").allocate(blitz::shape(3,2), {"three", "two"});

    auto &a1(bundle.array("a1"));
    EXPECT_EQ(0, a1.lbound(0));
    EXPECT_EQ(1, a1.ubound(0));
    EXPECT_EQ(0, a1.lbound(1));
    EXPECT_EQ(2, a1.ubound(1));
    EXPECT_EQ(6, a1.size());

    auto &a2(bundle.array("a2"));
    EXPECT_EQ(0, a2.lbound(1));
    EXPECT_EQ(1, a2.ubound(1));
    EXPECT_EQ(0, a2.lbound(0));
    EXPECT_EQ(2, a2.ubound(0));
    EXPECT_EQ(6, a2.size());

}
// -----------------------------------------------------------
#endif

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
