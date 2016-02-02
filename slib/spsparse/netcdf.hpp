#ifndef SPSPARSE_NETCDF_HPP
#define SPSPARSE_NETCDF_HPP

#include <functional>
#include <ibmisc/netcdf.hpp>
#include <spsparse/array.hpp>

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
	std::array<size_t, ArrayT::rank> shape;		// Extent of each dimension

	netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
	netCDF::NcVar vals_v = nc->getVar(vname + ".vals");

	std::vector<size_t> startp = {0, 0};		// SIZE, RANK
	std::vector<size_t> countp = {1, A->rank};	// Write RANK elements at a time
	for (auto ii = A->begin(); ii != A->end(); ++ii, ++startp[0]) {
		typename ArrayT::indices_type index = ii.index();
		typename ArrayT::val_type val = ii.val();

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
	std::array<size_t, AccumulatorT::rank> shape;		// Extent of each dimension

	netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
	netCDF::NcVar vals_v = nc->getVar(vname + ".vals");

	size_t size = vals_v.getDim(0).getSize();	// # non-zero elements

	std::vector<size_t> startp = {0, 0};		// SIZE, RANK
	std::vector<size_t> countp = {1, AccumulatorT::rank};	// Write RANK elements at a time
	std::array<typename AccumulatorT::index_type, AccumulatorT::rank> index;
	typename AccumulatorT::val_type val;

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
	std::vector<netCDF::NcDim> dims;		// Dimensions in NetCDF

	std::vector<size_t> dim_sizes;			// Length of our two dimensions.

	// Allocate the output, if we're reading
	if (ncio.rw == 'w') {
		dims = ibmisc::get_or_add_dims(ncio, dim_names, {A.size(), A.rank});

		auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});
		info_v.putAtt("shape", netCDF::ncUint64, A.rank, &A.shape[0]);

		get_or_add_var(ncio, vname + ".indices", netCDF::ncInt64, dims);
		get_or_add_var(ncio, vname + ".vals", netCDF::ncDouble, {dims[0]});
		ncio += std::bind(&nc_write_spsparse<ArrayT>, ncio.nc, &A, vname);
	} else {
		dims = ibmisc::get_dims(ncio, dim_names);

		// Read
		netCDF::NcVar info_v = ncio.nc->getVar(vname + ".info");
		auto shape_a = info_v.getAtt("shape");

		// Check the rank in NetCDF matches SpSparse rank
		size_t rank = shape_a.getAttLength();
		if (rank != ArrayT::rank) {
			(*spsparse_error)(-1,
				"Trying to read NetCDF sparse array of rank %ld into SpSparse array of rank %d",
				rank, ArrayT::rank);
		}

		if (alloc) {
			// Allocate + Read

			// Check the shape of the sparse array
			std::array<size_t, ArrayT::rank> shape;
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

}	// Namespace
#endif	// Guard
