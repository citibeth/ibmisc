// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <spsparse/xiter.hpp>
#include <iostream>

using namespace spsparse;

// The fixture for testing class Foo.
class XiterTest : public ::testing::Test {
protected:

	// You can do set-up work for each test here.
	XiterTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~XiterTest() {}

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


// Tests that STLXiter replicates standard
// STL iterator with new interface
TEST_F(XiterTest, STLXiter) {
	std::vector<int> vec = {0,2,4,6};
	auto ii(STLXiter<std::vector<int>::iterator>(vec.begin(), vec.end()));

	std::vector<int> new_vec;

	for (
		auto ii(STLXiter<std::vector<int>::iterator>(vec.begin(), vec.end()));
		!ii.eof(); ++ii) {
		new_vec.push_back(*ii);
	}
	EXPECT_EQ(vec, new_vec);
}

// Tests Join2Xiter
TEST_F(XiterTest, Join2Xiter) {
	std::vector<int> vec1 = {0,2,4,6};
	std::vector<int> vec2 = {0,1,2,3,4,5,6,7};
	std::vector<int> out;

	typedef STLXiter<std::vector<int>::iterator> XiterT;

	// -------- Test 1
	out.clear();
	for (auto ii(Join2Xiter<XiterT, XiterT>(
		XiterT(vec1.begin(), vec1.end()),
		XiterT(vec2.begin(), vec2.end())));
		!ii.eof(); ++ii)
	{
		EXPECT_EQ(*ii.i1, *ii.i2);
		out.push_back(*ii.i1);
	}
	EXPECT_EQ(out, std::vector<int>({0,2,4,6}));


	// -------- Test 2: Reverse order
	vec2 = {0,2,4,6};
	vec1 = {0,1,2,3,4,5,6,7};
	out.clear();
	for (auto ii(Join2Xiter<XiterT, XiterT>(
		XiterT(vec1.begin(), vec1.end()),
		XiterT(vec2.begin(), vec2.end())));
		!ii.eof(); ++ii)
	{
		out.push_back(*ii.i1);
	}
	EXPECT_EQ(out, std::vector<int>({0,2,4,6}));

	// -------- Test 3: More complex relations
	vec1 = {0,2,4,5,6,7,8,9};
	vec2 = {1,2,3,4,6};
	out.clear();
	for (auto ii(Join2Xiter<XiterT, XiterT>(
		XiterT(vec1.begin(), vec1.end()),
		XiterT(vec2.begin(), vec2.end())));
		!ii.eof(); ++ii)
	{
		out.push_back(*ii.i1);
	}
	EXPECT_EQ(out, std::vector<int>({2,4,6}));

}


// Tests Join3Xiter
TEST_F(XiterTest, Join3Xiter) {
	std::vector<int> vec1 = {0,2,4,6};
	std::vector<int> vec2 = {0,1,2,3,4,5,6,7};
	std::vector<int> vec3 = {1,2,3,6};
	std::vector<int> out;


	typedef STLXiter<std::vector<int>::iterator> XiterT;

	// -------- Test 1
	out.clear();
	for (auto ii(Join3Xiter<XiterT, XiterT, XiterT>(
		XiterT(vec1.begin(), vec1.end()),
		XiterT(vec2.begin(), vec2.end()),
		XiterT(vec3.begin(), vec3.end())));
		!ii.eof(); ++ii)
	{
		EXPECT_EQ(*ii.i1, *ii.i2);
		EXPECT_EQ(*ii.i1, *ii.i3);
		out.push_back(*ii.i1);
	}
	EXPECT_EQ(out, std::vector<int>({2,6}));

}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
