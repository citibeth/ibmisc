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
#include <ibmisc/filesystem.hpp>
#include <iostream>
#include <memory>

using namespace ibmisc;

// The fixture for testing class Foo.
class FilesystemTest : public ibmisc::Test {
protected:

    // You can do set-up work for each test here.
    FilesystemTest() : ibmisc::Test(false) {}    // keep=true
};

TEST_F(FilesystemTest, filesystem_simple)
{
    auto search_path(EnvSearchPath("PATH"));
    EXPECT_EQ("/bin/ls", search_path.locate("ls"));
}

// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
