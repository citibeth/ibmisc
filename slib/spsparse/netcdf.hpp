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

#ifndef SPSPARSE_NETCDF_HPP
#define SPSPARSE_NETCDF_HPP

#include <functional>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/netcdf.hpp>

namespace spsparse {

/** @defgroup netcdf netcdf.hpp
@brief Simple way to read/write SpSparse arrays via NetCDF.

@{
*/

template<class ArrayT>
void nc_write_spsparse(
    netCDF::NcGroup *nc,
    ArrayT *A,
    std::string const &vname);


template<class ArrayT>
void nc_write_spsparse(
    netCDF::NcGroup *nc,
    ArrayT *A,
    std::string const &vname)
{
    typedef std::array<typename ArrayT::index_type, ArrayT::rank> indices_type;

    std::array<size_t, ArrayT::rank> shape;     // Extent of each dimension

    netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
    netCDF::NcVar vals_v = nc->getVar(vname + ".vals");

    std::vector<size_t> startp = {0, 0};        // SIZE, RANK
    std::vector<size_t> countp = {1, A->rank};  // Write RANK elements at a time
    for (auto ii = A->begin(); ii != A->end(); ++ii, ++startp[0]) {
        auto index(ii->index());
        auto val(ii->value());

        indices_v.putVar(startp, countp, &index[0]);
        vals_v.putVar(startp, countp, &val);
    }
}
// --------------------------------------------------------

template<class AccumulatorT>
void nc_read_spsparse(
    netCDF::NcGroup *nc,
    AccumulatorT *A,
    std::string const &vname);

template<class AccumulatorT>
void nc_read_spsparse(
    netCDF::NcGroup *nc,
    AccumulatorT *A,
    std::string const &vname)
{
    std::array<size_t, AccumulatorT::rank> shape;       // Extent of each dimension

    netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
    netCDF::NcVar vals_v = nc->getVar(vname + ".vals");

    size_t size = vals_v.getDim(0).getSize();   // # non-zero elements

    std::vector<size_t> startp = {0, 0};        // SIZE, RANK
    std::vector<size_t> countp = {1, AccumulatorT::rank};   // Write RANK elements at a time
    std::array<typename AccumulatorT::index_type, AccumulatorT::rank> index;
    typename AccumulatorT::value_type val;

    for (; startp[0]<size; ++startp[0]) {
        indices_v.getVar(startp, countp, &index[0]);
        vals_v.getVar(startp, countp, &val);

        A->add(index, val);
    }
}

// -----------------------------------------------------------------
template<class ArrayT>
void ncio_spsparse(
    ibmisc::NcIO &ncio,
    ArrayT &A,
    bool alloc,
    std::string const &vname);

template<class ArrayT>
void ncio_spsparse(
    ibmisc::NcIO &ncio,
    ArrayT &A,
    bool alloc,
    std::string const &vname)
{
    std::vector<std::string> const dim_names({vname + ".size", vname + ".rank"});
    std::vector<netCDF::NcDim> dims;        // Dimensions in NetCDF

    std::vector<size_t> dim_sizes;          // Length of our two dimensions.

    // Allocate the output, if we're reading
    if (ncio.rw == 'w') {
        dims = ibmisc::get_or_add_dims(ncio, dim_names, {A.size(), A.rank});

        auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
        ibmisc::get_or_put_att(info_v, 'w', "shape", "uint64", A.shape);
//        info_v.putAtt("shape", ibmisc::nc_type("uint64"), A.rank, &A.shape[0]);

        get_or_add_var(ncio, vname + ".indices", "int64", dims);
        get_or_add_var(ncio, vname + ".vals", "double", {dims[0]});
        ncio += std::bind(&nc_write_spsparse<ArrayT>, ncio.nc, &A, vname);
    } else {
        dims = ibmisc::get_dims(ncio, dim_names);

        // Read
        netCDF::NcVar info_v = ncio.nc->getVar(vname + ".info");
        auto shape_a = info_v.getAtt("shape");

        // Check the rank in NetCDF matches SpSparse rank
        size_t rank = shape_a.getAttLength();
        if (rank != ArrayT::rank) {
            (*ibmisc::ibmisc_error)(-1,
                "Trying to read NetCDF sparse array of rank %ld into SpSparse array of rank %d",
                rank, ArrayT::rank);
        }

        if (alloc) {
            // Allocate + Read

            // Check the shape of the sparse array
            std::array<long, ArrayT::rank> shape;
            shape_a.getValues(&shape[0]);

            // Reserve space for the non-zero elements
            A.clear();
            A.set_shape(shape);
            A.reserve(ncio.nc->getDim(vname + ".size").getSize());
        }

        ncio += std::bind(&nc_read_spsparse<ArrayT>, ncio.nc, &A, vname);
    }
}



/** @} */

}   // Namespace
#endif  // Guard
