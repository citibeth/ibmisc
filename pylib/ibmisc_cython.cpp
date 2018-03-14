#include <cstdio>
#include <algorithm>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/cython.hpp>
#include <spsparse/SparseSet.hpp>
#include "ibmisc_cython.hpp"

static double const NaN = std::numeric_limits<double>::quiet_NaN();

namespace ibmisc {
namespace cython {


PyObject *linear_Weighted_shape(
    linear::Weighted &self)
{
    return py_build_tuple<pytype_long, long, 2>(self.shape());
}

PyObject *linear_Weighted_apply_weight(
    linear::Weighted &self,
    int dim,
    PyObject *A_s_py)            // A_b{nj_s} One row per variable
{
    auto shape(self.shape());

    // |j_s| = size of sparse input vector space (A_s)
    // |j_d] = size of dense input vector space (A_d)
    // |n| = number of variables being processed together

    // Convert input arrays
    auto A_s(np_to_blitz<double,2>(A_s_py, "A", {-1,-1}));    // Sparse indexed dense vector
    int n_n = A_s.extent(0);

    // Allocate output
    PyObject *B_s_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,1>{n_n});
    auto B_s(np_to_blitz<double,1>(B_s_py, "B", {-1}));
    B_s = NaN;

    // Run it!
    self.apply_weight(dim, A_s, B_s, true);

    return B_s_py;
}

PyObject *linear_Weighted_apply_M(
    linear::Weighted &self,
    PyObject *A_s_py,            // A_b{nj_s} One row per variable
    double fill,
    // std::string const &saccum_type,
    bool force_conservation)
{
    // Convert input arrays
    auto A_s(np_to_blitz<double,2>(A_s_py, "A", {-1,-1}));    // Sparse indexed dense vector
    int n_n = A_s.extent(0);

    // Convert accum_type
    // auto accum_type(parse_enum<AccumType>(saccum_type));

    // Allocate output
    auto shape(self.shape());
    PyObject *B_s_py = ibmisc::cython::new_pyarray<double,2>(
        std::array<int,2>{n_n, shape[0]});
    auto B_s(np_to_blitz<double,2>(B_s_py, "B_s_py", {-1,-1}));
    B_s = fill;

    // Run it!
    self.apply_M(A_s, B_s, linear::AccumType::REPLACE, force_conservation);

    return B_s_py;
}





}}    // namespace

