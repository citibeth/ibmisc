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
#ifndef IBMISC_PYLIB_EXAMPLES_HPP
#define IBMISC_PYLIB_EXAMPLES_HPP

namespace ibmisc {
namespace cython {

extern void cyexample_double_blitz(PyObject *a);

extern PyObject *cyexample_sparse_matrix();

extern std::unique_ptr<linear::Weighted> example_linear_weighted(
    std::string slinear_type);

}}    // namespace
#endif    // guard


