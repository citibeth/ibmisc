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

#include <iostream>
#include <ibmisc/VarTransformer.hpp>

namespace ibmisc {

const std::string VarTransformer::UNIT = "1";


void VarTransformer::set_dims(
    std::vector<std::string> const &outputs,
    std::vector<std::string> const &inputs,
    std::vector<std::string> const &scalars)
{
    std::array<std::vector<std::string> const *, RANK> source_dims {&outputs, &inputs, &scalars};

    // Copy the dimension labels, and add a "unit" dimension to the end
    for (int idim=0; idim<RANK; ++idim) {
        for (auto ii=source_dims[idim]->begin(); ii != source_dims[idim]->end(); ++ii)
            _dimensions[idim].insert(*ii);
        _dimensions[idim].insert(UNIT);
    }

    // Allocate the tensor to our size now.
    blitz::TinyVector<int, RANK> extent;
    for (int i=0; i<RANK; ++i) extent[i] = _dimensions[i].size();
    _tensor.reference(blitz::Array<double, RANK>(extent));
    _tensor = 0;
}



bool VarTransformer::set(std::string const &output, std::string const &input, std::string const &scalar, double val)
{
    bool is_good = true;
    if (!dim(OUTPUTS).contains(output)) {
        fprintf(stderr, "ERROR: VarTransformer::set(): output variable %s not defined.\n", output.c_str());
        is_good = false;
    }
    if (!dim(INPUTS).contains(input)) {
        fprintf(stderr, "ERROR: VarTransformer::set(): input variable %s not defined.\n", input.c_str());
        is_good = false;
    }
    if (!dim(SCALARS).contains(scalar)) {
        fprintf(stderr, "ERROR: VarTransformer::set(): scalar variable %s not defined.\n", scalar.c_str());
        is_good = false;
    }

    if (is_good) {
        int ioutput = dim(OUTPUTS).at(output);
        int iinput = dim(INPUTS).at(input);
        int iscalar = dim(SCALARS).at(scalar);
    
        _tensor(ioutput, iinput, iscalar) = val;
    }
    return is_good;
}

/** Local function for debugging output, not for re-use. */
static std::string dims_str(IndexSet<std::string> const &de)
{
    std::stringstream buf;
    int sz = de.size();
    for (int i=0; i<sz; ++i) buf << de[i] << ", ";
    return buf.str();
}

/** Instantiates the scalars with specific values, and returns a 2nd-order
matrix derived from the 3d-order tensor, in CSR format. */
CSRAndUnits VarTransformer::apply_scalars(
        std::vector<std::pair<std::string, double>> const &nvpairs)
{
#if 0
    printf("BEGIN VarTransformer::apply_scalars()\n");
    for (auto ele : nvpairs) {
        printf("    %s = %g\n", ele.first.c_str(), ele.second);
    }
#endif

//std::cout << "OUTPUTS = " << dims_str(dim(OUTPUTS)) << std::endl;
//std::cout << "INPUTS = " << dims_str(dim(INPUTS)) << std::endl;
//std::cout << "SCALARS = " << dims_str(dim(SCALARS)) << std::endl;


    int n_outputs_nu = dim(OUTPUTS).size()-1;       // # OUTPUTS no unit
    int n_inputs_wu = dim(INPUTS).size();
    int n_scalars_wu = dim(SCALARS).size(); // # SCALARS w/unit

    int unit_inputs = dim(INPUTS).size()-1;

    // Convert name/value pairs to a regular vector
    blitz::Array<double,1> scalars(n_scalars_wu);
    scalars = 0;
    for (auto ii = nvpairs.begin(); ii != nvpairs.end(); ++ii) {
        std::string const &nv_name = ii->first;
        double const val = ii->second;

        // If a provided scalar is not used for this VarTransformer, just ignore it.
        if (dim(SCALARS).contains(nv_name))
            scalars(dim(SCALARS).at(nv_name)) = val;
    }
    scalars(dim(SCALARS).at(UNIT)) = 1.0;

    // Take inner product of tensor with our scalars.
    CSRAndUnits ret(n_outputs_nu);
    for (int i=0; i < n_outputs_nu; ++i) {
        for (int j=0; j < n_inputs_wu; ++j) {
            double coeff = 0;
            for (int k=0; k < n_scalars_wu; ++k) {
                coeff += _tensor(i,j,k) * scalars(k);
            }

            // Output format: sparse matrix plus dense unit column
            if (j == unit_inputs) {
                ret.units[i] = coeff;
            } else {
                if (coeff != 0) ret.mat.add(i, j, coeff);
            }
        }
    }

//  std::cout << "apply_scalars() returning " << ret;

//  printf("END VarTransformer::apply_scalars()\n");
    return ret;
}

std::ostream &operator<<(std::ostream &out, CSRMatrix const &mat)
{
    out << "CSRMatrix :" << std::endl;
    for (unsigned int i=0; i < mat.size(); ++i) {
        out << "    " << i << ":";
        for (auto ele : mat[i]) {
            out << " (" << ele.second << "*[" << ele.first << "])";
        }
        out << std::endl;
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, CSRAndUnits const &matu)
{
    CSRMatrix const &mat(matu.mat);

    out << "CSRMatrix :" << std::endl;
    for (unsigned int i=0; i < mat.size(); ++i) {
        out << "    " << i << ": " << matu.units[i] << " +";
        for (auto ele : mat[i]) {
            out << " (" << ele.first << ", " << ele.second << ")";
        }
        out << std::endl;
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, VarTransformer const &vt)
{
    size_t n_outputs_nu = vt.dim(VarTransformer::OUTPUTS).size()-1;     // # OUTPUTS no unit
    size_t n_inputs_wu = vt.dim(VarTransformer::INPUTS).size();
    size_t n_scalars_wu = vt.dim(VarTransformer::SCALARS).size();   // # SCALARS w/unit

    size_t unit_outputs = vt.dim(VarTransformer::OUTPUTS).size()-1;
    size_t unit_inputs = vt.dim(VarTransformer::INPUTS).size()-1;
    size_t unit_scalars = vt.dim(VarTransformer::SCALARS).size()-1;


    for (int i=0; i<n_outputs_nu; ++i) {
        out << "    " << vt.dim(VarTransformer::OUTPUTS)[i] << " = ";

        // Count number of INPUTs used for this OUTPUT
        int nj = 0;
        std::vector<int> nk(n_inputs_wu, 0);
        for (int j=0; j < n_inputs_wu; ++j) {
            for (int k=0; k < n_scalars_wu; ++k) {
                double tijk = vt._tensor(i,j,k);
                if (tijk != 0) ++nk[j];
            }
            if (nk[j] > 0) ++nj;
        }

        // No RHS for this OUTPUT, quit
        if (nj == 0) {
            out << "0" << std::endl;
            continue;
        }

        // We DO have something on the RHS
        int jj = 0;
        for (int j=0; j < n_inputs_wu; ++j) {
            int nkj = nk[j];
            if (nkj == 0) continue;

            if (nkj > 1) out << "(";
            int kk = 0;
            for (int k=0; k < n_scalars_wu; ++k) {
                double val = vt._tensor(i,j,k);
                if (val == 0.0) continue;
                if (val != 1.0) out << val;
                if (k != unit_scalars) out << " " << vt.dim(VarTransformer::SCALARS)[k];

                if (kk != nkj-1) out << " + ";

                // Increment count of SEEN k values
                ++kk;
            }
            if (nkj > 1) out << ")";
            if (j != unit_inputs) out << " " << vt.dim(VarTransformer::INPUTS)[j];

            if (jj != nj-1) out << " + ";

//          out << vt.dim(VarTransformer::INPUTS)[j];

            // Increment count of SEEN j values
            ++jj;
        }
        out << std::endl;
    }
    return out;
}

}   // namespace ibmisc
