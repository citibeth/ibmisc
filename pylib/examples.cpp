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
