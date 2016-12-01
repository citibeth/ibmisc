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

#pragma once

#include <vector>
#include <array>
#include <blitz/array.h>
#include <ibmisc/IndexSet.hpp>

namespace ibmisc {

/** This is not officially part of the SparseMatrix.hpp framework.
Keeping it here is much simpler, even if not as general. */
class CSRMatrix : public std::vector<std::vector<std::pair<int, double>>> {
    typedef std::vector<std::vector<std::pair<int, double>>> super;
public:

    CSRMatrix(int nrow) : super(nrow) {}

    void add(int row, int col, double val)
        { (*this)[row].push_back(std::make_pair(col, val)); }

    friend std::ostream &operator<<(std::ostream &out, CSRMatrix const &mat);
};

/** A relation such as y = Ax + b can be represented as:
y = A'x' where x' = [x 1].  HOWEVER... such a representation
does not allow for easy application, unless a bunch of 1's are appended to x.
Therefore, it makes sense to decompose it to remove the "unit" element.
By convention, the unit element is the last one in a dimension.

This class stores the (smaller) vector x, along with the unit vector b. */
class CSRAndUnits {
public:
    CSRMatrix mat;
    std::vector<double> units;

    CSRAndUnits(int nrow) : mat(nrow), units(nrow, 0.0) {}

    friend std::ostream &operator<<(std::ostream &out, CSRAndUnits const &matu);
};

// -------------------------------------------------
// -------------------------------------------------

/** Inputs received from the GCM may not be in the same units, etc. as
the inputs to be sent to the ice model.  This class creates a matrix
to translate between those two.

We wish to make transformations such as the following:
   y1 = x1 * by_dt               (where by_dt is 1/ the timestep)
   y2 = x2 + 273.15
   y3 = x3 + i4

These transformations can be represented via the matrix:
   y = M [x 1]
where M is:
    x1      x2  x3  x4  unit
    -----------------------------
y1 |by_dt   0   0   0   0
y2 |0       1   0   0   273.15
y3 |0       0   1   1   0

Note that the Matrix M involves variables that won't be bound until
later (eg, by_dt, where dt is the (variable) COUPLING timestep).
Therefore, the we compute M by multiplying a 3d-order tensor T with the
vector [by_dt 1]:
    M = T . [by_dt 1]

Once all variables in [by_dt 1] are bound to values (eg, on each
coupling timestep), the matrix M may be computed.  M may then be used
to compute output values (y) based on input values (x).


Note that y and x will typically be bundles of vectors (for example,
y1 is SMB, y2 is energy flux, etc).  This computation needs to be computed
for each element of the y's and x's.

The dimensions of the 3d-order tensor T are called OUTPUTS, INPUTS and
SCALARS.  This class assigns strings to each element of each
dimension, thereby allowing the user to set elements of the tensor by
string.  This is inefficient, but is only done at initialization time.
And it prevents errors. */

class VarTransformer {
public:
    static const std::string UNIT;

    // Dimensions of our tensor (RANK is sentinel)
    enum {OUTPUTS, INPUTS, SCALARS, RANK};

protected:
    blitz::Array<double, RANK> _tensor;

    /** Name of each element in each dimension.
    NOTE: This includes a "unit" variable tacked to the end of each dimension. */

    std::array<IndexSet<std::string>,RANK> _dimensions;
public:

    /** Define the name of each element in a dimension. */
    void set_dims(
        std::vector<std::string> const &outputs,
        std::vector<std::string> const &inputs,
        std::vector<std::string> const &scalars);


    // Accessor methods...
    IndexSet<std::string> const &dim(int idim) const
        { return _dimensions[idim]; }

public:

    /** Set an element of the tensor, using name-based indexing.
    This is equivalent to:
            output += val * input * scalar

    @return true if all OK, false on error. */
    bool set(std::string const &output, std::string const &input, std::string const &scalar, double val);

    /** Instantiates the scalars with specific values, and returns a 2nd-order
    matrix derived from the 3d-order tensor, in CSR format. */
    CSRAndUnits apply_scalars(
        std::vector<std::pair<std::string, double>> const &nvpairs = {});

    friend std::ostream &operator<<(std::ostream &out, VarTransformer const &vt);
};

std::ostream &operator<<(std::ostream &out, CSRMatrix const &mat);
std::ostream &operator<<(std::ostream &out, CSRAndUnits const &matu);

/** Print out the tensor as readable symbolic equations.
Used to check and debug. */
std::ostream &operator<<(std::ostream &out, VarTransformer const &vt);


}
