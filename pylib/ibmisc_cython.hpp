#ifndef IBMISC_PYLIB_EXAMPLES_HPP
#define IBMISC_PYLIB_EXAMPLES_HPP

#include <Python.h>
#include <ibmisc/linear/linear.hpp>

namespace ibmisc {
namespace cython {

extern PyObject *linear_Weighted_shape(
    linear::Weighted *M);

extern PyObject *linear_Weighted_apply_weight(
    linear::Weighted *self,
    int dim,
    PyObject *A_s_py);            // A_b{nj_s} One row per variable

extern PyObject *linear_Weighted_apply_M(
    linear::Weighted *self,
    PyObject *A_s_py,            // A_b{nj_s} One row per variable
    double fill,
    // std::string const &saccum_type,
    bool force_conservation);


}}    // namespace
#endif    // guard
