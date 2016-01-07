#include <string>
#include <netcdf>
#include <sstream>
#include <ibmisc/netcdf.hpp>


using namespace netCDF;

namespace ibmisc {

/** Gets or creates an unlimited size dimension */
netCDF::NcDim get_or_add_dim(NcIO &ncio, std::string const &dim_name)
{
	bool err = false;

	NcDim dim = ncio.nc->getDim(dim_name);
	if (dim.isNull()){
		// The dim does NOT exist!
		if (ncio.rw == 'r') {
			(*ibmisc_error)(-1,
				"Dimension %s(unlimited) needs to exist when reading", dim_name.c_str());
		} else {
			// We're in write mode; make this dimension.
			return ncio.nc->addDim(dim_name);
		}
	}

	// Make sure existing dimension is unlimited
	if (!dim.isUnlimited()) {
		(*ibmisc_error)(-1, 
			"Attempt in get_or_add_dim() to change size from %d to unlimited",
			dim.getSize());
	}

	return dim;
}

netCDF::NcDim get_or_add_dim(NcIO &ncio, std::string const &dim_name, size_t dim_size)
{
	NcDim dim = ncio.nc->getDim(dim_name);
	if (dim.isNull()){
		// The dim does NOT exist!
		if (ncio.rw == 'r') {
			(*ibmisc_error)(-1,
				"Dimension %s(%s) needs to exist when reading", dim_name.c_str(), dim_size);
		} else {
			// We're in write mode; make this dimension.
			return ncio.nc->addDim(dim_name, dim_size);
		}
	}

	if (dim.getSize() != dim_size) {
		(*ibmisc_error)(-1, 
			"Attempt in get_or_add_dim() to change size from %ld to %ld",
			dim.getSize(), dim_size);
	}

	return dim;
}

std::vector<netCDF::NcDim> get_dims(
	NcIO &ncio,
	std::vector<std::string> const &sdims)
{
	size_t RANK = sdims.size();
	std::vector<netCDF::NcDim> ret(RANK);
	for (int k=0; k<RANK; ++k) ret[k] = ncio.nc->getDim(sdims[k]);
	return ret;
}

std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	std::vector<std::tuple<std::string, size_t>> const &sdims)
{
	size_t RANK = sdims.size();
	std::vector<netCDF::NcDim> ret(RANK);
	for (int k=0; k<RANK; ++k) {
		std::string const &dim_name(std::get<0>(sdims[k]));
		size_t const dim_size(std::get<1>(sdims[k]));

		if (dim_size < 0) {
			ret[k] = get_or_add_dim(ncio, dim_name);
		} else {
			ret[k] = get_or_add_dim(ncio, dim_name, dim_size);
		}
	}
	return ret;
}


}	// Namespace

