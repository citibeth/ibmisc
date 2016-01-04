#include <string>
#include <netcdf>
#include <sstream>


using namespace netCDF;

namespace ibmisc {

/** Gets or creates an unlimited size dimension */
netCDF::NcDim getOrAddDim(netCDF::NcGroup &nc, std::string const &dim_name)
{
	NcDim dim;
	bool err = false;

	try {
		// Look for the dim already existing
		dim = nc.getDim(dim_name);
	} catch(exceptions::NcBadDim) {
		// No existing dim, must create it.
		return nc.addDim(dim_name);
	}

	// Make sure existing dimension is unlimited
	if (!dim.isUnlimited()) {
		std::stringstream msg;
		msg << "Attempt in NcGroup::getOrAddDim to change size from " <<
			dim.getSize() << " to unlimited";
		throw exceptions::NcBadDim(msg.str(), __FILE__, __LINE__);
	}

	return dim;
}

netCDF::NcDim getOrAddDim(netCDF::NcGroup &nc, std::string const &dim_name, size_t dim_size)
{
	NcDim dim;
	try {
		// Look for the dim already existing
		dim = nc.getDim(dim_name);
	} catch(exceptions::NcBadDim) {
		// Must create the new dim
		return nc.addDim(dim_name, dim_size);
	}

	if (dim.getSize() != dim_size) {
		std::stringstream msg;
		msg << "Attempt in NcGroup::getOrAddDim to change size from " <<
			dim.getSize() << " to " << dim_size;
		throw exceptions::NcBadDim(msg.str(), __FILE__, __LINE__);
	}

	return dim;
}



}	// Namespace

