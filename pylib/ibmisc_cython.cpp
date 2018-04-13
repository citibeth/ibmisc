#include <cstdio>
#include <algorithm>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/cython.hpp>
#include <spsparse/SparseSet.hpp>
#include "ibmisc_cython.hpp"

static double const NaN = std::numeric_limits<double>::quiet_NaN();

namespace ibmisc {
namespace cython {


PyObject *linear_Weighted_type(
    linear::Weighted &self)
{
    return Py_BuildValue("s", self.type.str());
}

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


PyObject *linear_Weighted_to_coo(linear::Weighted const &self)
{
    int const rank = 2;


    // Allocate Python arrays and Blitzify
    std::array<PyObject *,rank> indices_np;
    std::array<blitz::Array<int,1>, rank> indices_bl;
    std::array<int,1> dims {(int)self.nnz()};
    for (int k=0; k<rank; ++k) {
        indices_np[k] = PyArray_FromDims(
            1, &dims[0], np_type_num<int>());
        indices_bl[k].reference(np_to_blitz<int,1>(indices_np[k], "index", dims));
    }
    PyObject *values_np = PyArray_FromDims(
        1, &dims[0], np_type_num<double>());
    auto values_bl(np_to_blitz<double,1>(values_np, "values", dims));


    // Fill with data
    self.to_coo(indices_bl[0], indices_bl[1], values_bl);


    // Bind it up in a tuple
    PyObject *indices_t = PyTuple_New(rank);
    PyObject *shape_t = PyTuple_New(rank);
    auto shape(self.shape());
    for (int k=0; k<rank; ++k) {
        // blitz::Array<typename ArrayT::val_type, rank> idx(indices(k));
        PyTuple_SetItem(indices_t, k, indices_np[k]);
        PyTuple_SetItem(shape_t, k, Py_BuildValue("i", shape[k]));
    }


    PyObject *data_t = Py_BuildValue("OO", values_np, indices_t);

    // For scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))
    PyObject *ret = Py_BuildValue("OO", data_t, shape_t);

    return ret;

}

PyObject *linear_Weighted_get_weights(linear::Weighted const &self, int idim)
{
    auto shape(self.shape());
    std::array<int,1> dims {(int)shape[0]};

    PyObject *W_np = PyArray_FromDims(
        1, &dims[0], np_type_num<double>());
    auto W(np_to_blitz<double,1>(W_np, "W", dims));

    self.get_weights(idim, W);
    return W_np;
}

}}    // namespace

