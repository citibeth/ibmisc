// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <spsparse/VectorCooArray.hpp>
#include <iostream>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif
#include <netcdf>
#include <spsparse/netcdf.hpp>

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

// The fixture for testing class Foo.
class SpSparseTest : public ::testing::Test {
protected:

	std::vector<std::string> tmpfiles;

	// You can do set-up work for each test here.
	SpSparseTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~SpSparseTest()
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


TEST_F(SpSparseTest, NetCDF) {
	// Make a simple VectorCooArray
	VectorCooArray<int, double, 2> arr1({5,6});
	arr1.add({1,2}, 2.);
	arr1.add({3,3}, 6.);
	arr1.add({4,5}, 1.);


	// If the constructor and destructor are not enough for setting up
	std::string fname("__netcdf_test.nc");
	tmpfiles.push_back(fname);
	::remove(fname.c_str());


	// Write
	{
		ibmisc::NcIO ncio(fname, NcFile::replace);
		ncio_spsparse(ncio, arr1, true, "arr1");
		ncio.close();
	}

	// Read with alloc
	VectorCooArray<int, double, 2> arr2;
	{
		ibmisc::NcIO ncio(fname, NcFile::read);
		ncio_spsparse(ncio, arr2, true, "arr1");
		ncio.close();
	}

	// Read without alloc
	VectorCooArray<int, double, 2> arr3(arr1.shape);
	{
		ibmisc::NcIO ncio(fname, NcFile::read);
		ncio_spsparse(ncio, arr3, false, "arr1");
		ncio.close();
	}


	// Check the three are equal
	EXPECT_EQ(arr1.size(), arr2.size());
	EXPECT_EQ(arr1.size(), arr3.size());
	for (auto ii1(arr1.begin()), ii2(arr2.begin()), ii3(arr3.begin()); ii1 != arr1.end(); ++ii1, ++ii2, ++ii3) {
		EXPECT_EQ(ii1.index(), ii2.index());
		EXPECT_EQ(ii1.val(), ii2.val());

		EXPECT_EQ(ii1.index(), ii3.index());
		EXPECT_EQ(ii1.val(), ii3.val());
	}

}


int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
	everytrace_init();
#endif
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
