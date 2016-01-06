#ifndef IBMISC_NCUTIL_HPP
#define IBMISC_NCUTIL_HPP

#include <netcdf>
#include <functional>

namespace ibmisc {

/** Creates a dimension if it doesn't already exist; or else checks
that the existing dimension has the requested size.
@param nc The NcGroup or NcFile to create the dimension in.
@param dim_name Name of dimension to create.
@param dim_size Size to create or check for.
@return The created/retrieved dimension.
*/
netCDF::NcDim getOrAddDim(netCDF::NcGroup &nc, std::string const &dim_name, size_t dim_size);

/** Creates an unlimited dimension if it doesn't already exist; or
else checks that the existing dimension has unlimited size.
@param nc The NcGroup or NcFile to create the dimension in.
@param dim_name Name of dimension to create.
@return The created/retrieved dimension.
*/
netCDF::NcDim getOrAddDim(netCDF::NcGroup &nc, std::string const &dim_name);

/** Used to keep track of future writes on NcDefine */
class NcWrites {
	std::vector<std::function<void ()>> _writes;
public:
	void operator+=(std::function<void ()> const &fn)
		{ _writes.push_back(fn); }

	void operator()()
	{
		for (auto ii=_writes.begin(); ii != _writes.end(); ++ii) (*ii)();
		_writes.clear();
	}
};

}	// Namespace
#endif	// Guard
