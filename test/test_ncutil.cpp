// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/ncutil.hpp>
#include <iostream>
#include <cstdio>
#include <netcdf>

using namespace ibmisc;
using namespace netCDF;

// The fixture for testing class Foo.
class NcutilTest : public ::testing::Test {
protected:

	std::vector<std::string> tmpfiles;

	// You can do set-up work for each test here.
	NcutilTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~NcutilTest()
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

//	  // The mock bar library shaed by all tests
//	  MockBar m_bar;
};

TEST_F(NcutilTest, getOrAddDim) {

	// If the constructor and destructor are not enough for setting up
	std::string fname("__ncutil_test.nc");
	tmpfiles.push_back(fname);
	::remove(fname.c_str());
	NcFile nc(fname, NcFile::replace);
	NcDim dim1 = getOrAddDim(nc, "dim1", 1);
	NcDim dim2 = getOrAddDim(nc, "dim2", 2);
	NcDim dimU = getOrAddDim(nc, "dimU");

	EXPECT_EQ(1, dim1.getSize());
	EXPECT_EQ(2, dim2.getSize());
	EXPECT_FALSE(dim1.isUnlimited());
	EXPECT_TRUE(dimU.isUnlimited());

	NcDim dim1b = getOrAddDim(nc, "dim1", 1);
	EXPECT_EQ(1, dim1b.getSize());
	EXPECT_FALSE(dim1b.isUnlimited());

	NcDim dimUb = getOrAddDim(nc, "dimU");
	EXPECT_TRUE(dimUb.isUnlimited());

	// Try to reset to unlimited
	EXPECT_THROW(getOrAddDim(nc, "dim1"), netCDF::exceptions::NcException);
	// Try to reset to another value
	EXPECT_THROW(getOrAddDim(nc, "dim1", 17), netCDF::exceptions::NcException);
	// Try to reset from unlimited
	EXPECT_THROW(getOrAddDim(nc, "dimU", 17), netCDF::exceptions::NcException);

}



int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
