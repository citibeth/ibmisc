#ifndef IBMISC_PYLIB_IBMISC_CYTHON_HPP
#define IBMISC_PYLIB_IBMISC_CYTHON_HPP

#include <Python.h>
#include <ibmisc/linear/linear.hpp>

namespace ibmisc {
namespace cython {

extern PyObject *linear_Weighted_type(
    linear::Weighted &self);

extern PyObject *linear_Weighted_shape(
    linear::Weighted &self);

extern PyObject *linear_Weighted_apply_weight(
    linear::Weighted &self,
    int dim,
    PyObject *A_s_py);            // A_b{nj_s} One row per variable

extern PyObject *linear_Weighted_apply_M(
    linear::Weighted &self,
    PyObject *A_s_py,            // A_b{nj_s} One row per variable
    double fill,
    // std::string const &saccum_type,
    bool force_conservation);

extern PyObject *linear_Weighted_to_coo(linear::Weighted const &self);

extern PyObject *linear_Weighted_get_weights(linear::Weighted const &self, int idim);

}}    // namespace
#endif    // guard
