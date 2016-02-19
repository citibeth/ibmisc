// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/iter.hpp>
#include <ibmisc/SegmentVector.hpp>
#include <iostream>
#include <cstdio>
#include <memory>
#include <map>

using namespace ibmisc;

// The fixture for testing class Foo.
class IterTest : public ::testing::Test {
protected:
	// You can do set-up work for each test here.
	IterTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~IterTest() {}

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


TEST_F(IterTest, deref_iter_unique_ptr)
{
	typedef std::vector<std::unique_ptr<int>> MyVector;
	typedef DerefRandomAccessIter<int, typename MyVector::iterator> iterator;
	typedef DerefRandomAccessIter<const int, typename MyVector::const_iterator> const_iterator;

	MyVector vec1;
	vec1.push_back(std::unique_ptr<int>(new int(1)));
	vec1.push_back(std::unique_ptr<int>(new int(2)));
	vec1.push_back(std::unique_ptr<int>(new int(3)));

	int i=1;
	auto i2 = iterator(vec1.begin());
	auto i1(iterator(vec1.begin()));

	for (auto ii=iterator(vec1.begin()); ii != iterator(vec1.end()); ++ii, ++i)
		EXPECT_EQ(*ii, i);

	i=1;
	for (const_iterator ii(vec1.cbegin()); ii != const_iterator(vec1.cend()); ++ii, ++i)
		EXPECT_EQ(*ii, i);

}
// -----------------------------------------------------------
TEST_F(IterTest, deref_iter_cptr)
{
	std::vector<int> ints(3);
	for (int i=0; i<3; ++i) ints[i] = i+1;

	typedef std::vector<int *> MyVector;
	typedef DerefRandomAccessIter<int, typename MyVector::iterator> iterator;
	typedef DerefRandomAccessIter<const int, typename MyVector::const_iterator> const_iterator;

	MyVector vec1;
	vec1.push_back(&ints[0]);
	vec1.push_back(&ints[1]);
	vec1.push_back(&ints[2]);

	int i=1;
	for (auto ii=iterator(vec1.begin()); ii != iterator(vec1.end()); ++ii, ++i)
		EXPECT_EQ(*ii, i);
}
// -----------------------------------------------------------
TEST_F(IterTest, second_iter)
{
	typedef std::map<int, double> MyMap;

	MyMap map1;
	map1.insert(std::make_pair(1, 1.0));
	map1.insert(std::make_pair(2, 2.0));
	map1.insert(std::make_pair(3, 3.0));

	typedef SecondIter<int, double, typename MyMap::iterator> iterator;
	typedef SecondIter<const int, const double, typename MyMap::const_iterator> const_iterator;

	int i=1;
	for (auto ii=iterator(map1.begin()); ii != iterator(map1.end()); ++ii, ++i) {
		EXPECT_EQ(*ii, (double)i);
		EXPECT_EQ(ii.key(), i);
	}
}
// -----------------------------------------------------------
TEST_F(IterTest, segment_vector)
{
	SegmentVector<std::string> svec(2);
	EXPECT_EQ(svec.end(), svec.begin());
	EXPECT_EQ(0, svec.size());


	svec.segment(1).push_back("unit");
	svec.segment(0).push_back("A");
	svec.segment(0).push_back("B");
	EXPECT_EQ(3, svec.size());

	EXPECT_EQ("A", svec[0]);
	EXPECT_EQ("B", svec[1]);
	EXPECT_EQ("unit", svec[2]);

	auto ii(svec.begin());
	EXPECT_EQ("A", *ii);
	++ii;
	EXPECT_EQ("B", *ii);
	++ii;
	EXPECT_EQ("unit", *ii);
	++ii;
	EXPECT_EQ(svec.end(), ii);
}
// -----------------------------------------------------------


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
