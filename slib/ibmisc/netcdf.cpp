#include <string>
#include <netcdf>
#include <sstream>
#include <ibmisc/netcdf.hpp>


using namespace netCDF;

namespace ibmisc {

void _check_nc_rank(
	netCDF::NcVar const &ncvar,
	int rank)
{
	// Check rank of NetCDF variable
	if (ncvar.getDimCount() != rank)
		(*ibmisc_error)(-1,
			"NetCDF variable of rank %ld does not match blitz::Array "
			"of rank %d", ncvar.getDimCount(), rank);
}

// ===============================================
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
	for (int k=0; k<RANK; ++k) {
		ret[k] = ncio.nc->getDim(sdims[k]);
		if (ret[k].isNull()) {
			(*ibmisc_error)(-1,
				"Dimension %s does not exist!", sdims[k].c_str());
		}
	}
	return ret;
}

std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	std::vector<std::string> const &dim_names,
	std::vector<size_t> const &dim_sizes)
{
	if (dim_names.size() != dim_sizes.size()) {
		(*ibmisc_error)(-1,
			"get_or_add_dims() requires dim_names[%ld] and dim_sizes[%ld] be of same length.",
			dim_names.size(), dim_sizes.size());
	}

	size_t RANK = dim_names.size();
	std::vector<netCDF::NcDim> ret(RANK);
	for (int k=0; k<RANK; ++k) {
		if (dim_sizes[k] < 0) {
			ret[k] = get_or_add_dim(ncio, dim_names[k]);
		} else {
			ret[k] = get_or_add_dim(ncio, dim_names[k], dim_sizes[k]);
		}
	}
	return ret;
}
// ====================================================
/** Error if var doesn't exist. */
netCDF::NcVar get_or_add_var(
	NcIO &ncio,
	std::string const &vname,
	netCDF::NcType const &nc_type,
	std::vector<netCDF::NcDim> const &dims)
{
	netCDF::NcVar ncvar;
	if (ncio.define) {
		ncvar = ncio.nc->addVar(vname, nc_type, dims);
	} else {
		ncvar = ncio.nc->getVar(vname);
		if (ncvar.isNull()) {
			(*ibmisc_error)(-1,
				"Variable %s required but not found", vname.c_str());
		}
	}
	return ncvar;
}

// ---------------------------------------------
/** Do linewrap for strings that are intended to be used as comment attributes in NetCDF files.
       see: http://www.cplusplus.com/forum/beginner/19034/
*/
std::string ncwrap( std::string const &str2, size_t width) {
	std::string str = "\n" + str2;
    size_t curWidth = width;
    while( curWidth < str.length() ) {
        std::string::size_type spacePos = str.rfind( ' ', curWidth );
        if( spacePos == std::string::npos )
            spacePos = str.find( ' ', curWidth );
        if( spacePos != std::string::npos ) {
            str[ spacePos ] = '\n';
            curWidth = spacePos + width + 1;
        }
    }

    return str;
}

}	// Namespace

