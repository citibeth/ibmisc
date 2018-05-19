#pragma once

#include <cstdio>
#include <gtest/gtest.h>
#include <map>

namespace ibmisc {

// The fixture for testing class Foo.
class Test : public ::testing::Test {
protected:

    bool keep;    // Remove files when we're done?
    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    Test(bool _keep) : keep(_keep) {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~Test()
    {
        if (!keep)
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

    /** Returns a temporary filename; will be deleted at the end. */
    std::string tmp_fname(std::string const &key)
    {
        tmpfiles.push_back(key);
        ::remove(key.c_str());
        return key;
    }


};

} // namespace
