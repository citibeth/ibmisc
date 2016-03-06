// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/memory.hpp>
#include <ibmisc/string.hpp>
#include <iostream>
#include <cstdio>
#include <memory>
#include <map>

using namespace std;
using namespace ibmisc;

// The fixture for testing class Foo.
class MemoryTest : public ::testing::Test {
protected:
    // You can do set-up work for each test here.
    MemoryTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~MemoryTest() {}

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

string hello = "hello";
string *get_hello()
{
    return &hello; 
}
unique_ptr<string> new_hello()
    { return unique_ptr<string>(new string("hello")); }

TEST_F(MemoryTest, lazy_ptr_test)
{
    LazyPtr<string> p0(get_hello());
    LazyPtr<string> p1(new_hello());
    LazyPtr<string> p2(&get_hello);
    LazyPtr<string> p3(&new_hello);

    EXPECT_EQ("hello", *p0);
    EXPECT_EQ("hello", *p1);
    EXPECT_EQ("hello", *p2);
    EXPECT_EQ("hello", *p3);
}

TEST_F(MemoryTest, string_printf)
{
    std::string s1(string_printf("Hello World"));
    EXPECT_EQ("Hello World", s1);

    std::string s2(string_printf("H%s%d", "ello ", 17));
    EXPECT_EQ("Hello 17", s2);
}


// -----------------------------------------------------------
// -----------------------------------------------------------


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
