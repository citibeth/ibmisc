#ifndef IBMISC_NETCDF_HPP
#define IBMISC_NETCDF_HPP

#include <netcdf>
#include <functional>
#include <tuple>
#include <memory>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/enum.hpp>

namespace ibmisc {

// ---------------------------------------------------
/** Used to keep track of future writes on NcDefine */
class NcIO {
	std::vector<std::function<void ()>> _io;
	netCDF::NcFile _mync;  // NcFile lacks proper move constructor
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
		_mync(filePath, fMode),
		own_nc(true),
		nc(&_mync),
		rw(fMode == netCDF::NcFile::FileMode::read ? 'r' : 'w'),
		define(rw == 'w') {}

	void operator+=(std::function<void ()> const &fn)
	{
		if (rw == 'r') fn();
		else _io.push_back(fn);
	}

	void operator()() {
		for (auto ii=_io.begin(); ii != _io.end(); ++ii) (*ii)();
		_io.clear();
	}

	void close() {
		if (own_nc) {
			(*this)();
			_mync.close();
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
	std::vector<std::string> const &dim_names);

extern std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	std::vector<std::string> const &dim_names,
	std::vector<size_t> const &dim_lens);

// ---------------------------------------------------------
template<class TypeT, int RANK>
std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	blitz::Array<TypeT, RANK> &val,
	std::vector<std::string> const &dim_names);

template<class TypeT, int RANK>
std::vector<netCDF::NcDim> get_or_add_dims(
	NcIO &ncio,
	blitz::Array<TypeT, RANK> &val,
	std::vector<std::string> const &dim_names)
{
	std::vector<size_t> dim_sizes(RANK);
	for (int k=0; k<RANK; ++k) dim_sizes[k] = val.extent(k);
	return get_or_add_dims(ncio, dim_names, dim_sizes);
}

// ===========================================================
// Variable Wrangling
netCDF::NcVar get_or_add_var(
	NcIO &ncio,
	std::string const &vname,
	netCDF::NcType const &nc_type,
	std::vector<netCDF::NcDim> const &dims);

template<class TypeT>
void get_or_put_var(netCDF::NcVar &ncvar, char rw,
	std::vector<size_t> const &startp,
	std::vector<size_t> const &countp,
	TypeT *dataValues);

template<class TypeT>
void get_or_put_var(netCDF::NcVar &ncvar, char rw,
	std::vector<size_t> const &startp,
	std::vector<size_t> const &countp,
	TypeT *dataValues)
{
	switch(rw) {
		case 'r' :
			ncvar.getVar(startp, countp, dataValues);
		break;
		case 'w' :
			ncvar.putVar(startp, countp, dataValues);
		break;
	}
}
// ========================================================
// Attribute Wrangling

template<class NcVarT, class AttrT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	const std::string &name, const netCDF::NcType &type,
	AttrT *data);

template<class NcVarT, class AttrT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	const std::string &name, const netCDF::NcType &type,
	AttrT *data, size_t len)
{
	switch(rw) {
		case 'w':
			ncvar.putAtt(name, type, len, data);
		break;
		case 'r':
			auto att(ncvar.getAtt(name));
			if (att.getAttLength() != len) {
				(*ibmisc_error)(-1,
					"Trying to read attribute %s of length %ld into C++ "
					"variable of length %ld",
					name.c_str(), att.getAttLength(), len);
			}
			att.getValues(data);
		break;
	}
}


template<class NcVarT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	std::string const &name,
	std::string &data);

template<class NcVarT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	std::string const &name,
	std::string &data)
{
	switch(rw) {
		case 'w':
			ncvar.putAtt(name, data);
		break;
		case 'r':
			auto att(ncvar.getAtt(name));
#if 0
			if (att.getAttLength() != len) {
				(*ibmisc_error)("Trying to read attribute %s of length %ld into C++ variable of length %ld", name, att.getAttLength(), len);
			}
#endif
			att.getValues(data);
		break;
	}
}

// ---------------------------------------

template<class NcVarT, class AttrT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	const std::string &name, const netCDF::NcType &type,
	AttrT &data)
{ get_or_put_att(ncvar, rw, name, type, &data, 1); }

// ---------------------------------------

template<class NcVarT, class AttrT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	const std::string &name, const netCDF::NcType &type,
	std::vector<AttrT> &data);

template<class NcVarT, class AttrT>
void get_or_put_att(
	NcVarT &ncvar, char rw,
	const std::string &name, const netCDF::NcType &type,
	std::vector<AttrT> &data)
{
	switch(rw) {
		case 'w':
			ncvar.putAtt(name, type, data.size(), &data[0]);
		break;
		case 'r':
			auto att(ncvar.getAtt(name));
			data.resize(att.getAttLength());
			att.getValues(&data[0]);
		break;
	}
}

// ---------------------------------------
// NcVarT = NcVar or NcGroup
// EnumT = boost::enum
template<class NcVarT, class EnumT>
void get_or_put_att_enum(
	NcVarT &ncvar, char rw,
	const std::string &name,
	EnumT &data)
{
	switch(rw) {
		case 'w':
			ncvar.putAtt(name, std::string(data.str()));
		break;
		case 'r':
			auto att(ncvar.getAtt(name));
			std::string sval;
			att.getValues(sval);
			data = parse_enum<EnumT>(sval);
		break;
	}
}
// ====================== Error Checking =============================
// ---------------------------------------------------
/** Check that blitz::Array is unit strides, column major.
For now, that's the only kind of Blitz variable we know how to write. */
template<class TypeT, int RANK>
static void _check_blitz_strides(blitz::Array<TypeT, RANK> const &val);

template<class TypeT, int RANK>
static void _check_blitz_strides(blitz::Array<TypeT, RANK> const &val)
{
	size_t expected_stride = 1;
	for (int k=RANK-1; k >= 0; --k) {
		if (val.stride(k) != expected_stride)
			(*ibmisc_error)(-1, "blitz::Array has unexpected stride, cannot read/write with NetCDF (for now).");
		expected_stride *= val.extent(k);
	}
}
// ---------------------------------------------------
/** Check that the NetCDF variable has the correct rank */
void _check_nc_rank(netCDF::NcVar const &ncvar, int rank);
// ---------------------------------------------------
template<class TypeT, int RANK>
static void _check_blitz_dims(
	netCDF::NcVar const &ncvar,
	blitz::Array<TypeT, RANK> const &val,
	char rw);

template<class TypeT, int RANK>
static void _check_blitz_dims(
	netCDF::NcVar const &ncvar,
	blitz::Array<TypeT, RANK> const &val,
	char rw)
{
	// Check dimensions of NetCDF var vs. blitz::Array
	for (int k=0; k<RANK; ++k) {
		netCDF::NcDim ncdim(ncvar.getDim(k));

		if ((rw == 'r' || !ncdim.isUnlimited()) &&
			(ncvar.getDim(k).getSize() != val.extent(k)))
		{
			(*ibmisc_error)(-1,
			"Dimension #%d (%ld) of blitz::Array must match %s:%s (%ld) in NetCDF",
			k, val.extent(k), ncvar.getName().c_str(),
			ncvar.getDim(k).getName().c_str(), ncvar.getDim(k).getSize());
		}
	}
}
// ---------------------------------------------------
template<class TypeT>
static void _check_vector_dims(
	netCDF::NcVar const &ncvar,
	std::vector<TypeT> const &val,
	char rw);

template<class TypeT>
static void _check_vector_size(
	netCDF::NcVar const &ncvar,
	std::vector<TypeT> const &val,
	char rw)
{
	// Check dimensions of NetCDF var vs. std::vector
	netCDF::NcDim ncdim(ncvar.getDim(0));

	if (rw == 'w' && ncdim.isUnlimited()) return;


	if (ncdim.getSize() != val.size()) {
		(*ibmisc_error)(-1,
			"Size (%ld) of std::vector must match %s:%s (%ld) in NetCDF",
			val.size(), ncvar.getName().c_str(),
			ncdim.getName().c_str(), ncdim.getSize());
	}
}
// ---------------------------------------------------
template<class TypeT, int RANK>
void nc_rw_blitz(
	netCDF::NcGroup *nc,
	char rw,
	blitz::Array<TypeT, RANK> *val,
	bool alloc,
	std::string const &vname)
{
	netCDF::NcVar ncvar = nc->getVar(vname);

	if (alloc && rw == 'r') {
		blitz::TinyVector<int,RANK> shape;
		for (int k=0; k<RANK; ++k) shape[k] = ncvar.getDim(k).getSize();
		val->resize(shape);
	}

	_check_blitz_strides(*val);
	_check_nc_rank(ncvar, RANK);
	_check_blitz_dims(ncvar, *val, rw);

	std::vector<size_t> startp(RANK);
	std::vector<size_t> countp(RANK);
	for (int k=0; k<RANK; ++k) {
		startp[k] = 0;	// Start on disk, which always starts at 0
		countp[k] = val->extent(k);
	}
	switch(rw) {
		case 'r' :
			ncvar.getVar(startp, countp, val->data());
		break;
		case 'w' :
			ncvar.putVar(startp, countp, val->data());
		break;
	}
}
// ---------------------------------------------------

template<class TypeT, int RANK>
blitz::Array<TypeT, RANK> nc_read_blitz(
	netCDF::NcGroup *nc,
	std::string const &vname)
{
	blitz::Array<TypeT, RANK> val;
	nc_rw_blitz(nc, 'r', &val, true, vname);
	return val;
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
	netCDF::NcVar ncvar = get_or_add_var(ncio, vname, nc_type, dims);

	// const_cast allows us to re-use nc_rw_blitz for read and write
	ncio += std::bind(&nc_rw_blitz<TypeT, RANK>,
		ncio.nc, ncio.rw, &val, alloc, vname);

}
// ----------------------------------------------------
// =================================================
// Specializations for std::vector instead of blitz::Array

/** Like get_or_add_dims() above, but specialized for std::vector */
template<class TypeT>
std::vector<netCDF::NcDim> get_or_add_dims(
	netCDF::NcGroup *nc,
	std::vector<TypeT> &val,
	std::array<std::string, 1> const &dim_names)
{
	std::vector<size_t> dim_sizes(1);
	dim_sizes[0] = val.size();
	return get_or_add_dims(nc, val, dim_names, dim_sizes);
}



template<class TypeT>
void nc_rw_vector(
	netCDF::NcGroup *nc,
	char rw,
	blitz::vector<TypeT> *val,
	bool alloc,
	std::string const &vname)
{

	netCDF::NcVar ncvar = nc->getVar(vname);

	if (alloc & rw == 'r') val->resize(ncvar.getDim(0).getSize());

	blitz::Array<TypeT, 1> bval(vector_blitz(*val));
	_check_nc_rank(ncvar, 1);
	_check_vector_dims(ncvar, *val, rw);

	std::vector<size_t> startp = {0};
	std::vector<size_t> countp = {val->size()};
	switch(rw) {
		case 'r' :
			ncvar.getVar(startp, countp, &val[0]);
		break;
		case 'w' :
			ncvar.putVar(startp, countp, &val[0]);
		break;
	}
}

template<class TypeT>
std::vector<TypeT> nc_read_vector(
	netCDF::NcGroup *nc,
	std::string const &vname)
{
	std::vector<TypeT> val;
	nc_rw_vector(nc, 'r', &val, true, vname);
	return val;
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
	get_or_add_var(ncio, vname, nc_type, dims);

	ncio += std::bind(&nc_rw_vector<TypeT>, ncio.nc, ncio.rw, &val, alloc, vname);
}
// ----------------------------------------------------

/** Do linewrap for strings that are intended to be used as comment attributes in NetCDF files.
       see: http://www.cplusplus.com/forum/beginner/19034/
*/
extern std::string ncwrap( std::string const &str, size_t width = 55 );





}	// Namespace
#endif	// Guard
