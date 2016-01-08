// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/netcdf.hpp>
#include <iostream>
#include <cstdio>
#include <netcdf>

using namespace ibmisc;
using namespace netCDF;

// The fixture for testing class Foo.
class NetcdfTest : public ::testing::Test {
protected:

	std::vector<std::string> tmpfiles;

	// You can do set-up work for each test here.
	NetcdfTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~NetcdfTest()
	{
		for (auto ii(tmpfiles.begin()); ii != tmpfiles.end(); ++ii) {
//			::remove(ii->c_str());
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

TEST_F(NetcdfTest, get_or_add_dim) {

	// If the constructor and destructor are not enough for setting up
	std::string fname("__netcdf_get_or_add_dim_test.nc");
	tmpfiles.push_back(fname);
	::remove(fname.c_str());
	NcIO ncio(fname, NcFile::replace);
	NcDim dim1 = get_or_add_dim(ncio, "dim1", 1);
	NcDim dim2 = get_or_add_dim(ncio, "dim2", 2);
	NcDim dimU = get_or_add_dim(ncio, "dimU");

	EXPECT_EQ(1, dim1.getSize());
	EXPECT_EQ(2, dim2.getSize());
	EXPECT_FALSE(dim1.isUnlimited());
	EXPECT_TRUE(dimU.isUnlimited());

	NcDim dim1b = get_or_add_dim(ncio, "dim1", 1);
	EXPECT_EQ(1, dim1b.getSize());
	EXPECT_FALSE(dim1b.isUnlimited());

	NcDim dimUb = get_or_add_dim(ncio, "dimU");
	EXPECT_TRUE(dimUb.isUnlimited());

	// Try to reset to unlimited
	EXPECT_THROW(get_or_add_dim(ncio, "dim1"), ibmisc::Exception);
	// Try to reset to another value
	EXPECT_THROW(get_or_add_dim(ncio, "dim1", 17), ibmisc::Exception);
	// Try to reset from unlimited
	EXPECT_THROW(get_or_add_dim(ncio, "dimU", 17), ibmisc::Exception);

}

TEST_F(NetcdfTest, blitz)
{
	std::string fname("__netcdf_blitz_test.nc");
	tmpfiles.push_back(fname);

	::remove(fname.c_str());
	
	blitz::Array<double,2> A(4,5);
	for (int i=0; i<A.extent(0); ++i) {
	for (int j=0; j<A.extent(1); ++j) {
	  A(i,j) = i*j;
	}}

	blitz::Array<double,2> B(A*2);

	// ---------- Write
	printf("Writing\n");
	{
	ibmisc::NcIO ncio(fname, NcFile::replace);
	auto dims = ibmisc::get_or_add_dims(ncio, A, {"dim4", "dim5"});
	ibmisc::ncio_blitz(ncio, A, true, "A", netCDF::ncDouble, dims);
	ibmisc::ncio_blitz(ncio, B, true, "B", netCDF::ncDouble, dims);
	ncio.close();
	}

	// ---------- Read
	printf("Reading\n");
	ibmisc::NcIO ncio(fname, NcFile::read);

	blitz::Array<double,2> A2, B2;
	auto dims = ibmisc::get_or_add_dims(ncio, A, {"dim4", "dim5"});
	ibmisc::ncio_blitz(ncio, A2, true, "A", netCDF::ncDouble, dims);
	ibmisc::ncio_blitz(ncio, B2, true, "B", netCDF::ncDouble, dims);
	ncio.close();

	for (int i=0; i<A.extent(0); ++i) {
	for (int j=0; j<A.extent(1); ++j) {
		EXPECT_EQ(A(i,j), A2(i,j));
	}}
//	std::cout << "A2" << A2 << std::endl;




}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
