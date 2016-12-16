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
#include <ibmisc/datetime.hpp>
#include <iostream>
#include <memory>

using namespace ibmisc;

// The fixture for testing class Foo.
class DatetimeTest : public ibmisc::Test {
protected:

    // You can do set-up work for each test here.
    DatetimeTest() : ibmisc::Test(false) {}    // keep=true
};

TEST_F(DatetimeTest, date)
{
    for (long j=0; j < 1000000; ++j) {
        JDate jda1(j);
        auto da(cal365.to_date(jda1));
        auto jda2(cal365.to_jdate(da));
        EXPECT_EQ(jda2, jda1);
    }

    Date da(2015,4,3);
    std::cout << da << std::endl;
    Date da2(cal365.to_date(cal365.to_jdate(da)));
    EXPECT_EQ(da2, da);

    Time tm(10,0,0);
    std::cout << tm << std::endl;
    Datetime dt(da,tm);
    std::cout << dt << std::endl;

    // Test operator=()
    da2 = da;
    EXPECT_EQ(da2, da);
    da2 = std::move(da);
}

TEST_F(DatetimeTest, time)
{
    for (long j=0; j < 86400L*1000000L; j += 1000289) {
        JTime jtm1(j);
        auto tm(to_time(jtm1));
        auto jtm2(to_jtime(tm));
        EXPECT_EQ(jtm2, jtm1);
    }

    Time tm(15,30,15);
    Time t2;
    t2 = tm;
    t2 = std::move(tm);

}

TEST_F(DatetimeTest, time_unit)
{
    TimeUnit tu(&cal365, Datetime(2015,3,1, 23,0,0,0), TimeUnit::SECOND);

    EXPECT_EQ(tu.to_datetime(0),            Datetime(2015,3,1, 23,0,0,0));
    EXPECT_EQ(tu.to_datetime(86400),        Datetime(2015,3,2, 23,0,0,0));
    EXPECT_EQ(tu.to_datetime(86400 + 1800), Datetime(2015,3,2, 23,30,0,0));
    EXPECT_EQ(tu.to_datetime(86400 + 3600), Datetime(2015,3,3, 0,0,0,0));
    EXPECT_EQ(tu.to_datetime(86400. + 2.*3600. + 1801.5), Datetime(2015,3,3, 1,30,1,500000));

//    TimeUnit tu2;
//    tu2 = tu;
//    tu2 = std::move(tu);
}

// -----------------------------------------------------------

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
