// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <spsparse/VectorCooArray.hpp>
#include <iostream>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

using namespace spsparse;
using namespace ibmisc;

// The fixture for testing class Foo.
class SpSparseTest : public ::testing::Test {
protected:

	// You can do set-up work for each test here.
	SpSparseTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~SpSparseTest() {}

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


TEST_F(SpSparseTest, VectorCooArray) {
	// Make a simple VectorCooArray
	VectorCooArray<int, double, 1> arr1({4});
	arr1.add({1}, 2.);
	arr1.add({3}, 6.);
	EXPECT_EQ(2, arr1.size());
	EXPECT_EQ(1, arr1.index(0, 0));
	EXPECT_EQ(3, arr1.index(0,1));
	EXPECT_EQ(2., arr1.val(0));

	// Test bounds checking
	try {
		arr1.add({17}, 4.);
		FAIL() << "Excpected spsparse::Exception";
	} catch(spsparse::Exception const &err) {
	} catch(...) {
		FAIL() << "Excpected spsparse::Exception";
	}

	// Test Move Constructor
	VectorCooArray<int, double, 1> arr2(std::move(arr1));
	EXPECT_EQ(2, arr2.size());
	EXPECT_EQ(1, arr2.index(0, 0));
	EXPECT_EQ(3, arr2.index(0,1));
	EXPECT_EQ(2., arr2.val(0));
}

/** Check that we can get sorted permutations properly. */
TEST_F(SpSparseTest, permutation) {
	VectorCooArray<int, double, 2> arr2({2,4});
	arr2.add({1,3}, 5.);
	arr2.add({1,2}, 3.);
	arr2.add({0,3}, 17.);

	std::vector<size_t> perm0(sorted_permutation(arr2, {0,1}));
	EXPECT_EQ(perm0, std::vector<size_t>({2,1,0}));

	std::vector<size_t> perm1(sorted_permutation(arr2, {1,0}));
	EXPECT_EQ(perm1, std::vector<size_t>({1,2,0}));

}


TEST_F(SpSparseTest, iterators) {
	VectorCooArray<int, double, 2> arr2({2,4});
	arr2.add({1,3}, 5.);
	arr2.add({1,2}, 3.);
	arr2.add({0,3}, 17.);
	arr2.add({1,2}, 15.);

	// Try out iterators a bit
	auto ii(arr2.begin());
	EXPECT_NE(ii, arr2.end());
	EXPECT_EQ(1, ii.index(0));
	EXPECT_EQ(3, ii.index(1));
	EXPECT_EQ(5., ii.val());
	ii.val() = 17.;
	EXPECT_EQ(17., ii.val());

	std::array<int,2> arr({1,3});
	EXPECT_EQ(arr, *ii);
	++ii;
	EXPECT_NE(ii, arr2.end());
	EXPECT_EQ(1, ii.index(0));
	EXPECT_EQ(2, ii.index(1));
	EXPECT_EQ(3., ii.val());

}

TEST_F(SpSparseTest, transpose) {
	// 2-D VectorCooArray; test consolidate
	VectorCooArray<int, double, 2> arr2({2,4});
	arr2.add({1,3}, 5.);
	arr2.add({1,2}, 3.);
	arr2.add({0,3}, 17.);
	arr2.add({0,1}, 14.);
	arr2.add({1,2}, 15.);

	arr2.transpose({0,1});
	EXPECT_EQ(std::vector<int>({1,1,0,0,1}), to_vector(arr2.indices(0)));
	EXPECT_EQ(std::vector<int>({3,2,3,1,2}), to_vector(arr2.indices(1)));	// j
	EXPECT_EQ(std::vector<double>({5., 3., 17., 14., 15.}), to_vector(arr2.vals()));

	arr2.transpose({1,0});
	EXPECT_EQ(std::vector<int>({3,2,3,1,2}), to_vector(arr2.indices(0)));	// j
	EXPECT_EQ(std::vector<int>({1,1,0,0,1}), to_vector(arr2.indices(1)));
	EXPECT_EQ(std::vector<double>({5., 3., 17., 14., 15.}), to_vector(arr2.vals()));

	arr2.transpose({1,0});
	EXPECT_EQ(std::vector<int>({1,1,0,0,1}), to_vector(arr2.indices(0)));
	EXPECT_EQ(std::vector<int>({3,2,3,1,2}), to_vector(arr2.indices(1)));	// j
	EXPECT_EQ(std::vector<double>({5., 3., 17., 14., 15.}), to_vector(arr2.vals()));

}


TEST_F(SpSparseTest, consolidate) {
	// 2-D VectorCooArray; test consolidate
	VectorCooArray<int, double, 2> arr2({2,4});
	arr2.add({1,3}, 5.);
	arr2.add({1,2}, 3.);
	arr2.add({0,3}, 17.);
	arr2.add({0,1}, 14.);
	arr2.add({1,2}, 15.);

	// Consolidate row major
	VectorCooArray<int, double, 2> arr3(arr2.shape);
	consolidate(arr3, arr2, {0,1});
	EXPECT_EQ(4, arr3.size());


	EXPECT_EQ(std::vector<int>({0,0,1,1}), to_vector(arr3.indices(0)));
	EXPECT_EQ(std::vector<int>({1,3,2,3}), to_vector(arr3.indices(1)));	// j
	EXPECT_EQ(std::vector<double>({14., 17., 18., 5.}), to_vector(arr3.vals()));

	EXPECT_EQ(std::vector<size_t>({0,2,4}), dim_beginnings(arr3));


	// Clear it out
	arr3.clear();
	EXPECT_EQ(0, arr3.size());

	// Consolidate column-major
	consolidate(arr3, arr2, {1,0});
	EXPECT_EQ(std::vector<int>({0,1,0,1}), to_vector(arr3.indices(0)));
	EXPECT_EQ(std::vector<int>({1,2,3,3}), to_vector(arr3.indices(1)));	// j
	EXPECT_EQ(std::vector<double>({14., 18., 17., 5.}), to_vector(arr3.vals()));
	EXPECT_EQ(std::vector<size_t>({0,1,2,4}), dim_beginnings(arr3));

}

TEST_F(SpSparseTest, dim_beginnings_iterators)
{
	typedef VectorCooArray<int, double, 2> VectorCooArrayT;
	VectorCooArrayT arr2({20,10});


	arr2.add({1,0}, 15.);
	arr2.add({1,3}, 17.);
	arr2.add({2,4}, 17.);
	arr2.add({6,4}, 10.);

	VectorCooArrayT arr3(arr2.shape);
	consolidate(arr3, arr2, {0,1});

	auto abegin(dim_beginnings(arr3));
	auto dbi(DimBeginningsXiter<VectorCooArrayT>(&arr3, 0, 1, abegin.begin(), abegin.end()));

	// First row
	EXPECT_EQ(1, *dbi);
	auto ii1(dbi.sub_xiter());
	EXPECT_EQ(0, *ii1);
	EXPECT_EQ(15., ii1.val());
	EXPECT_EQ(false, ii1.eof());
	++ii1;
	EXPECT_EQ(3, *ii1);
	EXPECT_EQ(17., ii1.val());
	EXPECT_EQ(false, ii1.eof());
	++ii1;
	EXPECT_EQ(true, ii1.eof());

	// Next row
	++dbi;
	EXPECT_EQ(2, *dbi);
	auto ii2(dbi.sub_xiter());
	EXPECT_EQ(4, *ii2);
	EXPECT_EQ(false, ii2.eof());
	++ii2;
	EXPECT_EQ(true, ii2.eof());

	// Next row
	++dbi;
	EXPECT_EQ(6, *dbi);
	auto ii6(dbi.sub_xiter());
	EXPECT_EQ(4, *ii6);
	EXPECT_EQ(false, ii6.eof());
	++ii6;
	EXPECT_EQ(true, ii6.eof());

}


TEST_F(SpSparseTest, dense)
{
	typedef VectorCooArray<int, double, 2> VectorCooArrayT;
	VectorCooArrayT arr2({20,10});

	arr2.add({1,0}, 15.);
	arr2.add({1,3}, 17.);
	arr2.add({2,4}, 17.);
	arr2.add({6,4}, 10.);

	blitz::Array<double, 2> dense(arr2.to_dense());
	int i,j;
	double sum=0;
	for (int i=0; i<20; ++i) {
	for (int j=0; j<10; ++j) {
		sum += dense(i,j);
	}}
	EXPECT_EQ(sum, 59.);

	for (auto ii(arr2.begin()); ii != arr2.end(); ++ii)
		EXPECT_EQ(dense(ii.index(0), ii.index(1)), ii.val());

}


TEST_F(SpSparseTest, dense_to_blitz)
{
	typedef VectorCooArray<int, double, 2> VectorCooArrayT;

	blitz::Array<double,2> dense1(4,5);
	dense1(2,3) = 5.0;
	dense1(2,4) = 6.0;
	dense1(0,1) = 7.0;

	VectorCooArrayT sparse1({4,5});
	to_sparse(sparse1, dense1);
	blitz::Array<double,2> dense2(sparse1.to_dense());

	for (int i=0; i<dense1.extent(0); ++i) {
	for (int j=0; j<dense1.extent(0); ++j) {
		double const d1(dense1(i,j));
		double const d2(dense2(i,j));

		EXPECT_EQ(d1, d2);
	}}

}


int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
	everytrace_init();
#endif
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
