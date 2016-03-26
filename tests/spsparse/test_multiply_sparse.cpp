/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// https://github.com/google/googletest/blob/master/googletest/docs/Primer.md

#include <iostream>
#include <random>
#include <gtest/gtest.h>
#include <ibmisc/blitz.hpp>
#include <spsparse/VectorCooArray.hpp>
#include <spsparse/multiply_sparse.hpp>
#include <spsparse/eigen.hpp>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

using namespace ibmisc;
using namespace spsparse;

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

//    // The mock bar library shaed by all tests
//    MockBar m_bar;
};


#if 0
TEST_F(SpSparseTest, VectorCooArray)
{
    typedef VectorCooArray<int, double, 2> VectorCooArrayT;

    VectorCooMatrix<int, double> row({2,10});
    row.add({0,8}, 6.);
    row.add({0,4}, 4.);
    row.add({0,0}, 2.);
    row.add({0,3}, 3.);
    row.add({1,8}, 3.);

    VectorCooVector<int, double> scale({10});
    scale.add({0}, 2.);
//  scale.add({3}, 3.);
    scale.add({4}, 4.);
    scale.add({8}, 4.);

    VectorCooMatrix<int, double> col({10,1});
    col.add({0,0}, 2.);
    col.add({3,0}, 3.);
//  col.add({4,0}, 4.);
    col.add({8,0}, 5.);


    // ------------------ Do the same single row-col test, but with "full" matrix mutliplier
    VectorCooVector<int, double> eye({10});
    for (int i=0; i<10; ++i) eye.add({i}, 1.);

    VectorCooArray<int, double, 2> ret2;
    multiply(ret2, 1.0,
        &eye, row, '.', &scale, col, '.', &eye);
std::cout << "ret2: ";
std::cout << ret2 << std::endl;

    EXPECT_EQ(2, ret2.size());
    EXPECT_EQ(std::vector<int>({0,1}), to_vector(ret2.indices(0)));
    EXPECT_EQ(std::vector<int>({0,0}), to_vector(ret2.indices(1))); // j
    EXPECT_EQ(std::vector<double>({128., 60.}), to_vector(ret2.vals()));
}

#endif

// ------------------------------------------------------
void test_random_MM_multiply(unsigned int dsize, int seed)
{
    std::default_random_engine generator(seed);
    auto dim_distro(std::bind(std::uniform_int_distribution<int>(0,dsize-1), generator));
    auto val_distro(std::bind(std::uniform_real_distribution<double>(0,1), generator));

    VectorCooMatrix<int, double> A({dsize,dsize});
    VectorCooMatrix<int, double> B({dsize,dsize});
    int nranda = (int)(val_distro() * (double)(dsize*dsize));
    for (int i=0; i<nranda; ++i) A.add({dim_distro(), dim_distro()}, val_distro());
    int nrandb = (int)(val_distro() * (double)(dsize*dsize));
    for (int i=0; i<nrandb; ++i) B.add({dim_distro(), dim_distro()}, val_distro());

    VectorCooVector<int, double> eye({dsize});
    for (int k=0; k<dsize; ++k) eye.add({k}, 1.0);

// std::cout << "A: " << A << std::endl;
// std::cout << "B: " << B << std::endl;

    VectorCooMatrix<int, double> C;
    multiply(C,1.0,
        (VectorCooVector<int, double> *)0,  // scalei
        A, '.',
        &eye,   // scalej
//      (VectorCooVector<int, double> *)0,  // scalej
        B, '.',
        (VectorCooVector<int, double> *)0); // scalek

// std::cout << "C: " << C << std::endl;

    // --------- Compare to dense matrix multiplication
    auto Ad(A.to_dense());
    auto Bd(B.to_dense());
    auto Cd(C.to_dense());

    double usum = 0;
    for (int i=0; i<dsize; ++i) {
    for (int j=0; j<dsize; ++j) {
        double sum=0;
        for (int k=0; k<dsize; ++k) {
            sum += Ad(i,k) * Bd(k,j);
        }
        EXPECT_DOUBLE_EQ(sum, Cd(i,j));
        usum += sum;
    }}
    printf("MM seed = %d  sizes = [%ld, %ld, %ld]  usum = %f\n", seed, A.size(), B.size(), C.size(), usum);
}

TEST_F(SpSparseTest, random_MM_multiply)
{
    for (int seed=1; seed<1000; ++seed)
        test_random_MM_multiply(5,seed);
}
// ---------------------------------------------------------
void test_random_MV_multiply(unsigned int dsize, int seed)
{
    std::default_random_engine generator(seed);
    auto dim_distro(std::bind(std::uniform_int_distribution<int>(0,dsize-1), generator));
    auto val_distro(std::bind(std::uniform_real_distribution<double>(0,1), generator));

    VectorCooMatrix<int, double> A({dsize,dsize});
    VectorCooVector<int, double> B({dsize});
    int nranda = (int)(val_distro() * (double)(dsize*dsize));
    for (int i=0; i<nranda; ++i) A.add({dim_distro(), dim_distro()}, val_distro());
    int nrandb = (int)(val_distro() * (double)dsize);
    for (int i=0; i<nrandb; ++i) B.add({dim_distro()}, val_distro());

    VectorCooVector<int, double> eye({dsize});
    for (int k=0; k<dsize; ++k) eye.add({k}, 1.0);

    VectorCooVector<int, double> C;
    multiply(C,1.0,
        (VectorCooVector<int, double> *)0,  // scalei
        A, '.',
//      &eye,   // scalej
        (VectorCooVector<int, double> *)0,  // scalej
        B);

    // --------- Compare to dense matrix multiplication
    auto Ad(A.to_dense());
    auto Bd(B.to_dense());
    auto Cd(C.to_dense());

    double usum = 0;
    bool err = false;
    for (int i=0; i<dsize; ++i) {
        double sum=0;
        for (int k=0; k<dsize; ++k) {
            sum += Ad(i,k) * Bd(k);
        }
        if (sum != Cd(i)) {
            printf("(%d): %g vs %g\n", i,sum, Cd(i));

            FAIL() << "Error on multiply" << std::endl
                << "    A: " << A << std::endl
                << "    B: " << B << std::endl
                << "    C: " << C << std::endl;
            goto break_mult;
        }
        usum += sum;
    }
break_mult: ;
}

TEST_F(SpSparseTest, random_MV_multiply)
{
    for (int seed=1; seed<1000; ++seed) {
        test_random_MV_multiply(5,seed);
    }
}
// -----------------------------------------------------------
void test_random_VV_multiply(unsigned int dsize, int seed)
{
    std::default_random_engine generator(seed);
    auto dim_distro(std::bind(std::uniform_int_distribution<int>(0,dsize-1), generator));
    auto val_distro(std::bind(std::uniform_real_distribution<double>(0,1), generator));

    VectorCooVector<int, double> A({dsize});
    VectorCooVector<int, double> B({dsize});
    int nranda = (int)(val_distro() * (double)dsize);
    for (int i=0; i<nranda; ++i) A.add({dim_distro()}, val_distro());
    int nrandb = (int)(val_distro() * (double)dsize);
    for (int i=0; i<nrandb; ++i) B.add({dim_distro()}, val_distro());

    A.consolidate({0});
    B.consolidate({0});

    VectorCooVector<int, double> C({dsize});
    multiply_ele(C, A, B);

    // --------- Compare to dense matrix multiplication
    auto Ad(A.to_dense());
    auto Bd(B.to_dense());
    auto Cd(C.to_dense());

    double usum = 0;
    bool err = false;
    for (int i=0; i<dsize; ++i) {
        double prod = Ad(i) * Bd(i);
        if (Cd(i) != prod) {
            printf("(%d): %g vs %g\n", i, prod, Cd(i));

            FAIL() << "Error on multiply" << std::endl
                << "    A: " << A << std::endl
                << "    B: " << B << std::endl
                << "    C: " << C << std::endl;
            goto break_mult;
        }
        usum += prod;
    }
break_mult:
    printf("MV seed = %d  sizes = [%ld, %ld, %ld]  usum = %f\n", seed, A.size(), B.size(), C.size(), usum);
}

TEST_F(SpSparseTest, random_VV_multiply)
{
    for (int seed=1; seed<1000; ++seed) {
        test_random_VV_multiply(5,seed);
    }
}
// -----------------------------------------------------------
TEST_F(SpSparseTest, eigen_conv)
{
    size_t dsize = 17;
    SparseSet<int,int> dimi;
    SparseSet<int,int> dimj;
    SparseTriplets<VectorCooMatrix<int, double>> A({&dimi, &dimj});
    A.set_shape({dsize, dsize});
    A.add({4,5}, 1.);
    A.add({4,1}, 2.);
    A.add({0,1}, 3.);

    auto Ae(A.to_eigen());

    VectorCooMatrix<int, double> A2;
    eigenM_to_sparseM(A2, Ae, {&dimi, &dimj});

    A.M.consolidate({0,1});
    A2.consolidate({0,1});

    auto Ai0(to_vector(A.M.indices(0)));
    auto Ai1(to_vector(A.M.indices(1)));
    auto Av(to_vector(A.M.vals()));
    EXPECT_EQ(Ai0, to_vector(A2.indices(0)));
    EXPECT_EQ(Ai1, to_vector(A2.indices(1)));
    EXPECT_EQ(Av, to_vector(A2.vals()));
}
// -----------------------------------------------------------
void test_random_MM_multiply_eigen(unsigned int dsize, int seed)
{
    std::default_random_engine generator(seed);
    auto dim_distro(std::bind(std::uniform_int_distribution<int>(0,dsize-1), generator));
    auto val_distro(std::bind(std::uniform_real_distribution<double>(0,1), generator));

    SparseSet<int,int> dimi;
    SparseSet<int,int> dimj;
    SparseSet<int,int> dimk;

    SparseTriplets<VectorCooMatrix<int, double>> A({&dimi, &dimj});
    A.set_shape({dsize,dsize});
    SparseTriplets<VectorCooMatrix<int, double>> B({&dimj, &dimk});
    B.set_shape({dsize,dsize});

    int nranda = (int)(val_distro() * (double)(dsize*dsize));
    for (int i=0; i<nranda; ++i) A.add({dim_distro(), dim_distro()}, val_distro());
    int nrandb = (int)(val_distro() * (double)(dsize*dsize));
    for (int i=0; i<nrandb; ++i) B.add({dim_distro(), dim_distro()}, val_distro());

    Eigen::SparseMatrix<double> Ae, Be, Ce;
    Ae = A.to_eigen('.');
    Be = B.to_eigen('.');
    Ce = Ae * Be;

    VectorCooMatrix<int, double> C({A.M.shape[0], B.M.shape[1]});
    eigenM_to_sparseM(C, Ce, {&dimi, &dimk});

    // --------- Compare to dense matrix multiplication
    auto Ad(A.M.to_dense());
    auto Bd(B.M.to_dense());
    auto Cd(C.to_dense());
    double usum = 0;
    for (int i=0; i<dsize; ++i) {
    for (int j=0; j<dsize; ++j) {
        double sum=0;
        for (int k=0; k<dsize; ++k) {
            sum += Ad(i,k) * Bd(k,j);
        }
        EXPECT_DOUBLE_EQ(sum, Cd(i,j));
        usum += sum;
    }}
//  printf("MM seed = %d  sizes = [%ld, %ld, %ld]  usum = %f\n", seed, A.M.size(), B.M.size(), C.size(), usum);
}

TEST_F(SpSparseTest, random_MM_multiply_eigen)
{
    for (int seed=1; seed<1000; ++seed)
        test_random_MM_multiply_eigen(5,seed);
}
// ---------------------------------------------------------
int main(int argc, char **argv) {
#ifdef USE_EVERYTRACE
    everytrace_init();
#endif
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
