#include <ibmisc/cython.hpp>
#include <spsparse/VectorCooArray.hpp>

#include "examples.hpp"

namespace ibmisc {
namespace cython {

static void example_double_blitz(blitz::Array<double,1> &a)
{
    for (int i=a.lbound(0); i<=a.ubound(0); ++i) {
        a(i) *= 2;
    }
}

/** Cython interface... converts np.array to blitz::Array */
void cyexample_double_blitz(PyObject *a)
{
    blitz::Array<double,1> ab(np_to_blitz<double,1>(a, "A", {-1}));
    example_double_blitz(ab);
}

PyObject *cyexample_sparse_matrix()
{
    spsparse::VectorCooArray<long, double, 2> M({4,6});
    M.add({2,3}, 2.2);
    M.add({0,5}, 1.1);
    PyObject *tuple = spsparse_to_tuple(M);
    return tuple;
}

}}
