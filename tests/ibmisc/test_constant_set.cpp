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
#include <ibmisc/ConstantSet.hpp>
#include <ibmisc/string.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include <map>

using namespace std;
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
}

#if 0
TEST_F(ConstantSetTestNcIO, constant_set)
{

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
#else
std::string const constants_cdl =
    "netcdf __constants2 {\n"
    "dimensions:\n"
    "\tone = 1 ;\n"
    "variables:\n"
    "\tdouble constant.lhe(one) ;\n"
    "\t\tconstant.lhe:value = 2500000. ;\n"
    "\t\tconstant.lhe:units = \"J kg-1\" ;\n"
    "\t\tconstant.lhe:description = \"latent heat of evap at 0 C 2.5008d6\" ;\n"
    "\tdouble constant.lhm(one) ;\n"
    "\t\tconstant.lhm:value = 334000. ;\n"
    "\t\tconstant.lhm:units = \"J kg-1\" ;\n"
    "\t\tconstant.lhm:description = \"latent heat of melt at 0 C 334590\" ;\n"
    "data:\n"
    "\n"
    " constant.lhe = _ ;\n"
    "\n"
    " constant.lhm = _ ;\n"
    "\n"
    "}\n";




TEST_F(ConstantSetTest, constant_set_ncread)
{
    // Test the NetCDF I/O
    std::string fname("__constants2");
//    tmpfiles.push_back(fname+".cdl");
//    tmpfiles.push_back(fname+".nc");
    ::remove((fname+".cdl").c_str());

    // Write .cdl to temporary file
    {ofstream fout;
        fout.open((fname+".cdl").c_str());
        fout << constants_cdl;
        fout.close();
    }

    // Convert to .nc
     std::string const cmd(string_printf("ncgen -o %s.nc %s.cdl", fname.c_str(), fname.c_str()));

    FILE *in = popen(cmd.c_str(), "r");
    EXPECT_NE(nullptr, in);
 
    char buf[512];
    while(fgets(buf, sizeof(buf), in)!=NULL){
        std::cout << buf;
    }
    pclose(in);


    // Read the .nc file
    UTSystem ut_system("");
    ConstantSet cs;
    cs.init(&ut_system);
    {NcIO ncio(fname+".nc", 'r');
        cs.read_nc(ncio.nc, "");
    }

    double const lhe = cs.get_as("constant.lhe", "J kg-1");
    EXPECT_EQ(2500000., lhe);
    EXPECT_EQ(334000., cs.get_as("constant.lhm", "J kg-1"));
}
#endif

// -----------------------------------------------------------


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

