#include <ibmisc/cython.hpp>

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

}}
