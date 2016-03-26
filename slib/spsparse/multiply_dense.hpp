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


/** @defgroup multiply_dense multiply_dense.hpp
@brief Multiplication of sparse matrices by dense vectors.

@{
*/

#if 0

/** (Sparse Matrix) * (Dense Vector) */
template<class VectorCooMatrixT, class AccumulatorT>
void multiply(
    VectorCooMatrixT const &M,
    blitz::Array<double,1> const &x,
    AccumulatorT &y,
    bool handle_nan = false,
    bool transpose = false)
{
    for (auto ii = M.begin(); ii != M.end(); ++ii) {
        std::array<IndexT,2> index(ii.index());
        if (transpose) std::swap(index[0], index[1]);
        double val = ii.val() * x(index[1]);
        if (!handle_nan || !(std::isnan(val) || std::isinf(val))) {
            y.add(  {{index[0]}},  val);
        }
    }
}

template<class VectorCooMatrixT, class AccumulatorT>
void multiplyT(
    VectorCooMatrixT const &M,
    blitz::Array<double,1> const &x,
    AccumulatorT &y,
    bool handle_nan = false)
{ return multiply(M, x, y, handle_nan, true); }
#endif

// ---------------------------------------------------

/** @} */
