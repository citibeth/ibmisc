/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <ibmisc/fortranio.hpp>
#include <spsparse/eigen.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <ibmisc/linear/compressed.hpp>

using namespace std;
using namespace ibmisc;
using namespace blitz;
using namespace spsparse;
using namespace std::placeholders;


// -----------------------------------------
// ============== From icebin/eigen_types.hpp
// Types that will be used throughout as template arguments
typedef long sparse_index_type;
typedef int dense_index_type;
typedef double val_type;

// -----------------------------------------
typedef spsparse::MakeDenseEigen<sparse_index_type, val_type, 0, dense_index_type> MakeDenseEigenT;
template<int RANK>
    using TupleListT = MakeDenseEigenT::TupleListT<RANK>;
template<int RANK>
    using DenseArrayT = blitz::Array<val_type,RANK>;
typedef MakeDenseEigenT::SparseSetT SparseSetT;
typedef MakeDenseEigenT::EigenSparseMatrixT EigenSparseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic> EigenRowVectorT;
typedef Eigen::DiagonalMatrix<val_type, Eigen::Dynamic> EigenDiagonalMatrixT;
// -----------------------------------------


extern "C" void write_hntr40(
    int const &ima, int const &jma, float const &offia, float const &dlata,
    int const &imb, int const &jmb, float const &offib, float const &dlatb,
    float const &datmis);

// The fixture for testing class Foo.
class LinearTest : public ::testing::Test {
protected:

    std::vector<std::string> tmpfiles;

    // You can do set-up work for each test here.
    LinearTest() {}

    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~LinearTest()
    {
        for (auto ii(tmpfiles.begin()); ii != tmpfiles.end(); ++ii) {
//          ::remove(ii->c_str());
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
};

void sample_makedense(MakeDenseEigenT::AccumT &&accum)
{

    accum.add({0,0},.5);
    accum.add({10,10},.5);
    accum.add({20,20},.5);


    //accum.add({30,30},.5);    // Nullspace
    accum.add({00,40},.5);
    accum.add({10,40},.5);
    accum.add({20,40},.5);
    //accum.add({30,40},.5);
}

/** Tests the overlap matrix (area of exchange grid) conforms to
some basic properties such a matrix should have. */
TEST_F(LinearTest, overlap)
{
    // ============ BvA1
    // Compute an Eigen-based overlap matrix
    std::array<MakeDenseEigenT::SparseSetT, 2> dims;

    dims[0].set_sparse_extent(40);
    dims[1].set_sparse_extent(50);

    // Add extra unneeded dimensions
    dims[0].add_dense(15);
    dims[1].add_dense(15);

    MakeDenseEigenT BvA_o(
        std::bind(&sample_makedense, _1),
        {SparsifyTransform::ADD_DENSE},
        {&dims[0], &dims[1]}, '.');

    // Add extra unneeded dimensions
    dims[0].add_dense(25);
    dims[1].add_dense(25);

 
    auto BvA(BvA_o.to_eigen());

    // Produce a scaled regridding matrix
    auto sBvA(sum(BvA, 0, '-'));

    // Scale and stick into a linear::Weighted_Eigen
    linear::Weighted_Eigen BvA1;
//    BvA1.tmp.take(BvA1.dims[0], std::move(dims[0]));
//    BvA1.tmp.take(BvA1.dims[1], std::move(dims[1]));
    BvA1.dims[0] = &dims[0];
    BvA1.dims[1] = &dims[1];
    BvA1.M.reset(new MakeDenseEigenT::EigenSparseMatrixT(map_eigen_diagonal(sBvA) * BvA));
    BvA1.wM.reference(sum(BvA,0,'+'));
    BvA1.Mw.reference(sum(BvA,1,'+'));

    // Make it not conservative
    BvA1.wM(2) *= .9;
    BvA1.conservative = false;

    // Store in NetCDF
    std::string fname("__linear1.nc");
    tmpfiles.push_back(fname);
    ::remove(fname.c_str());
    {NcIO ncio(fname, 'w');
        dims[0].ncio(ncio, "dimB");
        dims[1].ncio(ncio, "dimA");
        BvA1.ncio(ncio, "BvA");
    }

    // ============ BvA2: Eigen, read from disk
    std::array<MakeDenseEigenT::SparseSetT, 2> dims2;
    linear::Weighted_Eigen BvA2;
        BvA2.dims[0] = &dims2[0];
        BvA2.dims[1] = &dims2[1];
    {NcIO ncio(fname, 'r');
        // Read back NetCDF
        dims2[0].ncio(ncio, "dimB");
        dims2[1].ncio(ncio, "dimA");
        BvA2.ncio(ncio, "BvA");
    }

    // =========== BvA2x: Compressed
    linear::Weighted_Compressed BvA2x(compress(BvA2));

    // Store in NetCDF
    std::string fname3("__linear3.nc");
    tmpfiles.push_back(fname3);
    ::remove(fname3.c_str());
    {NcIO ncio(fname3, 'w');
        BvA2x.ncio(ncio, "BvA");
    }

    // =========== BvA3: Read it back
    linear::Weighted_Compressed BvA3;
    {NcIO ncio(fname3, 'r');
        BvA3.ncio(ncio, "BvA");
    }


    // =================== Compare Indices
    linear::Weighted const * const BvA1p = &BvA1;
    linear::Weighted const * const BvA3p = &BvA3;


#if 0
printf("--------- BvA3.M 1\n");
for (auto ii(BvA3.M.generator()); ++ii; ) {
    printf("    {%d %d} %g\n", ii->index(0), ii->index(1), ii->value());
}
#endif

    // ---------------- Matrix
    {
        auto ii1(begin(*BvA1.M));
        auto ii3(BvA3.M.generator());
        int nnz = 0;
        for (;;) {
            // While...
            bool const end1 = (ii1 == end(*BvA1.M));
            bool const end3 = (!++ii3);
            EXPECT_EQ(end1,end3);
            if (end1) break;
            if (end1 != end3) break;    // Don't torture the loop any more if something went wrong

            // Loop body
            EXPECT_EQ(BvA1.dims[0]->to_sparse(ii1->index(0)), ii3->index(0));
            EXPECT_EQ(BvA1.dims[1]->to_sparse(ii1->index(1)), ii3->index(1));
            EXPECT_EQ(ii1->value(), ii3->value());

            // Increment
            ++nnz;
            ++ii1;
        }
        EXPECT_EQ(BvA3.M.nnz(), nnz);
        EXPECT_EQ(6,nnz);
    }

    // ------------------ Weight Vectors
    for (int idim=0; idim<2; ++idim) {
        auto &weights1(idim == 0 ? BvA1.wM : BvA1.Mw);
        auto &weights3(BvA3.weights[idim]);

        int ii1 = 0;
        auto ii3(weights3.generator());
        int nnz = 0;
        for (;; ++nnz, ++ii1) {
            // Skip over "accidental" nullspace
            while (weights1(ii1) == 0 && ii1 < weights1.extent(0)) ++ii1;

            // While...
            bool const end1 = (ii1 == weights1.extent(0));
            bool const end3 = (!++ii3);
            EXPECT_EQ(end1,end3);
            if (end1) break;
            if (end1 != end3) break;    // Don't torture the loop any more if something went wrong

            // Loop body
            EXPECT_EQ(BvA1.dims[idim]->to_sparse(ii1), ii3->index(0));
            EXPECT_EQ(weights1(ii1), ii3->value());
        }
        EXPECT_EQ(BvA3.weights[idim].nnz(), nnz);
    }

    // -------- Misc.
    EXPECT_EQ(BvA3p->conservative, BvA1p->conservative);
    auto shape1(BvA1p->shape());
    auto shape3(BvA3p->shape());
    EXPECT_EQ(shape1[0], shape3[0]);
    EXPECT_EQ(shape1[1], shape3[1]);


    // ================== Compare Multiplication (without conservation correction)
    for (int force_conservation=0; force_conservation<2; ++force_conservation) {
        printf("============== force_conservation = %d\n", force_conservation);
        int const nk = 1;
        blitz::Array<double,2> aa(nk,shape1[1]);
        for (int i=0; i<aa.extent(1); ++i) aa(0,i) = 32*i*i - 17;

        blitz::Array<double,2> bb1(nk,shape1[0]);
        blitz::Array<double,2> bb3(nk,shape1[0]);

        bb1 = -17;
        bb3 = -17;

        BvA1p->apply_M(aa, bb1, linear::AccumType::REPLACE, force_conservation);
        BvA3p->apply_M(aa, bb3, linear::AccumType::REPLACE, force_conservation);

        for (int k=0; k<nk; ++k) {    // Per var
            for (int i=0; i<bb1.extent(1); ++i) {
                //printf("(%d, %d): (%g, %g)\n", k,i,bb1(k,i),bb3(k,i));
                EXPECT_DOUBLE_EQ(bb1(k,i), bb3(k,i));
            }
        }
    }

    // NOT TESTED:
    //    Other AccumTypes

}


int main(int argc, char **argv) {
    // For test, configure Everytrace to silently throw exceptions (which we can catch)
    everytrace_init();
//    everytrace_exit = &everytrace_exit_silent_exception;

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
