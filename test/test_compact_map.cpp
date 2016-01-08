// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/CompactMap.hpp>
#include <iostream>
#include <cstdio>

using namespace ibmisc;

// The fixture for testing class Foo.
class CompactMapTest : public ::testing::Test {
protected:

	std::vector<std::string> tmpfiles;

	// You can do set-up work for each test here.
	CompactMapTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~CompactMapTest()
	{
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

class Value {
	std::string const _name;
public:
	int const val;
	std::string const &name() const { return _name; }
	Value(std::string const &name, int _val) : _name(name), val(_val) {}

	// http://www.learncpp.com/cpp-tutorial/93-overloading-the-io-operators/
	friend std::ostream &operator<<(std::ostream &out, Value const &con);
};

std::ostream &operator<<(std::ostream &out, Value const &val)
{
	out << "<" << val.name() << "=" << val.val << ">";
	return out;
}




TEST_F(CompactMapTest, compact_map) {
	CompactMap<std::string, Value> cm;

	EXPECT_EQ(0, cm.insert(Value("val0", 0)));
	Value val1("val1", 1);
	EXPECT_EQ(1, cm.insert(val1));
	EXPECT_EQ(2, cm.insert(Value("val3", 3)));
	EXPECT_EQ(3, cm.insert(Value("val2", 2)));
	EXPECT_EQ(4, cm.size());

	std::cout << "Example of operator<<():" << std::endl;
	std::cout << cm << std::endl;

	EXPECT_EQ(2, cm.index("val3"));
	EXPECT_EQ("val3", cm["val3"].name());
	EXPECT_EQ("val3", cm[2].name());

	{
	auto ii(cm.begin());
	EXPECT_EQ("val0", ii->name());
	EXPECT_EQ(0, ii->val); ++ii;
	EXPECT_EQ(1, ii->val); ++ii;
	EXPECT_EQ(3, ii->val); ++ii;
	EXPECT_EQ(2, ii->val); ++ii;
	EXPECT_EQ(ii, cm.end());
	}

	// Const iterators
	{
	CompactMap<std::string, Value> const &ccm(cm);

	auto ii(ccm.begin());
	EXPECT_EQ("val0", ii->name());
	EXPECT_EQ(0, ii->val); ++ii;
	EXPECT_EQ(1, ii->val); ++ii;
	EXPECT_EQ(3, ii->val); ++ii;
	EXPECT_EQ(2, ii->val); ++ii;
	EXPECT_EQ(ii, cm.end());
	}


}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
