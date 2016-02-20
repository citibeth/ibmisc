// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <gtest/gtest.h>
#include <ibmisc/VarTransformer.hpp>
#include <iostream>
#include <cstdio>
#include <memory>
#include <map>

using namespace ibmisc;

// The fixture for testing class Foo.
class VarTransformerTest : public ::testing::Test {
protected:
	std::vector<std::string> tmpfiles;

	// You can do set-up work for each test here.
	VarTransformerTest() {}

	// You can do clean-up work that doesn't throw exceptions here.
	virtual ~VarTransformerTest()
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

TEST_F(VarTransformerTest, segment_vector)
{
	IndexSet<std::string> dim_outputs({"len[cm]", "T[F]", "total_mass[kg]"});
	IndexSet<std::string> dim_inputs({"len[in]", "T[C]", "mass_per_timestep[kg s-1]"});
	IndexSet<std::string> dim_scalars({"dt[s]"});

	EXPECT_EQ(3, dim_outputs.size());
	EXPECT_EQ(3, dim_outputs.keys().size());

	VarTransformer vt;
	vt.set_dims(dim_outputs.keys(), dim_inputs.keys(), dim_scalars.keys());

	bool ok = true;
	ok = ok && vt.set("len[cm]", "len[in]", "1", 2.54);
	ok = ok && vt.set("T[F]", "T[C]", "1", 9./5.);
	ok = ok && vt.set("T[F]", "1", "1", 32.);
	ok = ok && vt.set("total_mass[kg]", "mass_per_timestep[kg s-1]", "dt[s]", 1.);
	EXPECT_TRUE(ok);

#if 0
	printf("------------------- vt =\n");
	std::cout << vt << std::endl;
	printf("-------------------\n");
#endif

	ibmisc::CSRAndUnits trans = vt.apply_scalars({std::make_pair("dt[s]", 17.0)});

#if 0
	printf("------------------- trans =\n");
	std::cout << trans << std::endl;
	printf("-------------------------\n");
#endif

	// Try a unit conversion with our newly minted matrix!
	blitz::Array<double,1> oval(dim_outputs.size());
	blitz::Array<double,1> ival(dim_inputs.size());
	ival(dim_inputs.at("len[in]")) = 12.;
	ival(dim_inputs.at("T[C]")) = 10.;
	ival(dim_inputs.at("mass_per_timestep[kg s-1]")) = 3.;

	// Try it out...
	for (int xi=0; xi<dim_outputs.size(); ++xi) {
		double sum = 0;
		std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
		for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
			int xj = xjj->first;
			double io_val = xjj->second;

			sum += io_val * ival(xj);
		}
		oval(xi) = sum + trans.units[xi];
	}

	EXPECT_DOUBLE_EQ(30.48, oval(0));
	EXPECT_DOUBLE_EQ(50., oval(1));
	EXPECT_DOUBLE_EQ(51., oval(2));

#if 0
	for (size_t i=0; i<dim_inputs.size(); ++i) {
		printf("ival[%s] = %g\n", dim_inputs[i].c_str(), ival(i));
	}
	for (size_t i=0; i<dim_outputs.size(); ++i) {
		printf("oval[%s] = %g\n", dim_outputs[i].c_str(), oval(i));
	}
#endif




}

// -----------------------------------------------------------


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
