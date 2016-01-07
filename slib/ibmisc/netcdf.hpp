#ifndef IBMISC_NCUTIL_HPP
#define IBMISC_NCUTIL_HPP

#include <netcdf>
#include <functional>
#include <memory>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/blitz.hpp>

namespace ibmisc {

// ---------------------------------------------------
/** Used to keep track of future writes on NcDefine */
class NcIO {
	std::vector<std::function<void ()>> _io;
	std::unique_ptr<netCDF::NcFile> _mync;  // NcFile lacks proper move constructor
	bool own_nc;
public:
	netCDF::NcGroup * const nc;
	char const rw;
	const bool define;

	/** @param _mode:
        'd' : Define and write (if user calls operator() later)
        'w' : Write only
		'r' : Read only (if user calls operator() later) */
	NcIO(netCDF::NcGroup *_nc, char _mode) :
		nc(_nc),
		own_nc(false),
		rw(_mode == 'd' ? 'w' : 'r'),
		define(_mode == 'd') {}

	NcIO(std::string const &filePath, netCDF::NcFile::FileMode fMode) :
		_mync(new netCDF::NcFile(filePath, fMode)),
		own_nc(true),
		nc(_mync.get()),
		rw(fMode == netCDF::NcFile::FileMode::read ? 'r' : 'w'),
		define(rw == 'w') {}

	void operator+=(std::function<void ()> const &fn)
		{ _io.push_back(fn); }

	void operator()() {
		for (auto ii=_io.begin(); ii != _io.end(); ++ii) (*ii)();
		_io.clear();
	}

	void close() {
		if (own_nc) {
			(*this)();
			_mync->close();
		} else {
			(*ibmisc_error)(-1, "NcIO::close() only valid on NcGroups it owns.");
		}
	}
};
// ===========================================================
// Dimension Wrangling
// ---------------------------------------------------
/** Creates a dimension if it doesn't already exist; or else checks
that the existing dimension has the requested size.
@param nc The NcGroup or NcFile to create the dimension in.
@param dim_name Name of dimension to create.
@param dim_size Size to create or check for.
@return The created/retrieved dimension.
*/
netCDF::NcDim get_or_add_dim(NcIO &ncio, std::string const &dim_name, size_t dim_size);

/** Creates an unlimited dimension if it doesn't already exist; or
else checks that the existing dimension has unlimited size.
@param nc The NcGroup or NcFile to create the dimension in.
@param dim_name Name of dimension to create.
@return The created/retrieved dimension.
*/
netCDF::NcDim get_or_add_dim(NcIO &ncio, std::string const &dim_name);

// ---------------------------------------------------
/** Convert dimensions from strings to NcDim */
extern std::vector<netCDF::NcDim> get_dims(
	NcIO &ncio,
	std::vector<std::string> const &sdims);

extern std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	std::vector<std::tuple<std::string, size_t>> const &ssdims);

// ---------------------------------------------------------
template<class TypeT, int RANK>
std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	blitz::Array<TypeT, RANK> &val,
	std::array<std::string, RANK> const &sdims);

template<class TypeT, int RANK>
std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	blitz::Array<TypeT, RANK> &val,
	std::array<std::string, RANK> const &sdims)
{
	std::vector<std::tuple<std::string, size_t>> ssdims(RANK);
	for (int k=0; k<RANK; ++k) {
		std::get<0>(ssdims[k]) = sdims[k];
		std::get<1>(ssdims[k]) = val.extent(k);
	}
	return get_or_add_dims(ncio, ssdims);
}

// ========================================================
// ---------------------------------------------------
template<class TypeT, int RANK>
static void _check_blitz(
	netCDF::NcVar const &ncvar,
	blitz::Array<TypeT, RANK> const &val,
	char rw)
{
	// Check that blitz::Array is unit strides, column major
	size_t expected_stride = 1;
	for (int k=RANK-1; k >= 0; --k) {
		if (val.stride(k) != expected_stride)
			(*ibmisc_error)(-1, "blitz::Array has unexpected stride, cannot read/write with NetCDF (for now).");
		expected_stride *= val.extent(k);
	}

	// Check rank of NetCDF variable
	if (ncvar.getDimCount() != RANK)
		(*ibmisc_error)(-1, "NetCDF variable of rank %ld does not match blitz::Array of rank %d", ncvar.getDimCount(), RANK);


	// Check dimensions of NetCDF var vs. blitz::Array
	for (int k=0; k<RANK; ++k) {
		netCDF::NcDim ncdim(ncvar.getDim(k));

		if (   (rw == 'r' || !ncdim.isUnlimited())  &&
			(ncvar.getDim(k).getSize() != val.extent(k))
		) {
			(*ibmisc_error)(-1,
			"Dimension #%d (%ld) of blitz::Array must match %s (%ld) in NetCDF",
			k, val.extent(k), ncvar.getDim(k).getName().c_str(), ncvar.getDim(k).getSize());
		}
	}

}


template<class TypeT, int RANK>
void nc_rw_blitz(
	netCDF::NcGroup *nc,
	char rw,
	blitz::Array<TypeT, RANK> *val,
	std::string const &vname)
{
	netCDF::NcVar ncvar = nc->getVar(vname);
	_check_blitz(ncvar, *val, rw);

	std::vector<size_t> startp(RANK);
	std::vector<size_t> endp(RANK);
	for (int k=0; k<RANK; ++k) {
		startp[k] = 0;	// Start on disk, which always starts at 0
		endp[k] = val->extent(k);
	}
	switch(rw) {
		case 'r' :
			ncvar.getVar(startp, endp, val->data());
		break;
		case 'w' :
			ncvar.putVar(startp, endp, val->data());
		break;
	}
}

/** Define and write a blitz::Array. */
template<class TypeT, int RANK>
void ncio_blitz(
	NcIO &ncio,
	blitz::Array<TypeT, RANK> &val,
	bool alloc,
	std::string const &vname,
	netCDF::NcType const &nc_type,
	std::vector<netCDF::NcDim> const &dims)
{
	if (alloc && ncio.rw == 'r') {
		blitz::TinyVector<int,RANK> shape;
		for (int k=0; k<RANK; ++k) shape[k] = dims[k].getSize();
		val.resize(shape);
	}

	netCDF::NcVar ncvar;
	if (ncio.define) {
		ncvar = ncio.nc->addVar(vname, nc_type, dims);
	} else {
		ncvar = ncio.nc->getVar(vname);
	}
	_check_blitz(ncvar, val, ncio.rw);

	// const_cast allows us to re-use nc_rw_blitz for read and write
	ncio += std::bind(&nc_rw_blitz<TypeT, RANK>,
		ncio.nc, ncio.rw, &val, vname);

}
// ----------------------------------------------------
// =================================================
// Specializations for std::vector instead of blitz::Array

/** Like get_or_add_dims() above, but specialized for 
template<class TypeT>
void std::vector<netCDF::NcDim> get_or_add_dims(
	netCDF::NcGroup *nc,
	std::vector<TypeT> &val,
	std::array<std::string, 1> const &sdims)
{
	std::vector<std::tuple<std::string, size_t>> ssdims(RANK);
	std::get<0>(ssdims[0]) = sdims[0];
	std::get<1>(ssdims[0]) = val.extent(k);
	return get_or_add_dims(nc, val, ssdims);
}

/** Define and write a std::vector. */
template<class TypeT>
void ncio_vector(
	NcIO &ncio,
	std::vector<TypeT> &val,
	bool alloc,			// Should we allocate val?
	std::string const &vname,
	netCDF::NcType const &nc_type,
	std::vector<netCDF::NcDim> const &dims)
{
	if (alloc & ncio.rw == 'r') val.resize(dims[0].getSize());
	ncio_blitz(ncio, vector_to_blitz(val), false, vname, nc_type, dims);
}
// ----------------------------------------------------






}	// Namespace
#endif	// Guard
