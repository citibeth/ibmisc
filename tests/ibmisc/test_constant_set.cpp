// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/ConstantSet.hpp>
#include <iostream>
#include <cstdio>
#include <memory>
#include <map>

using namespace ibmisc;
using namespace netCDF;

// The fixture for testing class Foo.
class ConstantSetTest : public ::testing::Test {
protected:
    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    ConstantSetTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~ConstantSetTest()
    {
        for (auto ii(tmpfiles.begin()); ii != tmpfiles.end(); ++ii) {
            ::remove(ii->c_str());
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

TEST_F(ConstantSetTest, units)
{
    UTSystem ut_system("");
    ut_system.parse("m");
}

TEST_F(ConstantSetTest, constant_set)
{
    UTSystem ut_system("");

    ConstantSet cs;
    cs.init(&ut_system);
    EXPECT_EQ(0, cs.size());

    cs.set("length", 1700., "cm", "Length of our thing");
    cs.set("width", 4., "m", "Width of our thing");
    EXPECT_EQ(2, cs.size());

    auto ii(cs.index.begin());
    EXPECT_EQ(ii.index(), 0);
    EXPECT_EQ(*ii, "length");
    ++ii;
    EXPECT_EQ(ii.index(), 1);
    EXPECT_EQ(*ii, "width");
    ++ii;
    EXPECT_EQ(cs.index.end(), ii);

    EXPECT_EQ("cm", cs[0].units);
    EXPECT_EQ("m", cs[1].units);

    double length_m = cs.get_as("length", "m");
    EXPECT_DOUBLE_EQ(17., length_m);
    double width_m = cs.get_as("width", "m");
    EXPECT_DOUBLE_EQ(4., width_m);

// This was never implemented...
//  // Try out operator<<
//  printf("ConstantSet: ");
//  std::cout << cs << std::endl;

    // Test the NetCDF I/O
    std::string fname("__constants.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());

    {NcIO ncio(fname, NcFile::replace);
    cs.ncio(ncio, "constants");}

    ConstantSet cs2;
    cs2.init(&ut_system);
    {NcIO ncio(fname, NcFile::read);
    cs2.ncio(ncio, "constants");}

    EXPECT_EQ(cs.size(), cs2.size());
    for (auto ii(cs.index.begin()); ii != cs.index.end(); ++ii) {
        size_t ix = ii.index();
        auto datum(cs[ix]);
        auto datum2(cs2[ix]);

        EXPECT_EQ(datum.name, datum2.name);
        EXPECT_EQ(datum.units, datum2.units);
        EXPECT_EQ(datum.description, datum2.description);
        EXPECT_EQ(datum.val, datum2.val);
    }
}

// -----------------------------------------------------------


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
