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

#ifndef IBMISC_NETCDF_HPP
#define IBMISC_NETCDF_HPP

#include <typeinfo>
#include <netcdf>
#include <boost/any.hpp>
#include <functional>
#include <tuple>
#include <memory>
#include <blitz/array.h>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/enum.hpp>
#include <ibmisc/memory.hpp>
#include <type_traits>

namespace ibmisc {

extern bool netcdf_debug;

// ---------------------------------------------------
// Convert template types to NetCDF types

// http://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf.html
// (byte, char, short, ushort, int, uint, int64, uint64, float or real, double, string) 
template<class T>
inline std::string get_nc_type()
{
    (*ibmisc_error)(-1,
        "get_nc_type(): Unknown type: '%s'", typeid(T).name());
}

template<> inline std::string get_nc_type<bool>()
    { return "char"; }

template<> inline std::string get_nc_type<int8_t>()
    { return "byte"; }

template<> inline std::string get_nc_type<char>()
    { return "char"; }

template<> inline std::string get_nc_type<int16_t>()
    { return "short"; }

template<> inline std::string get_nc_type<uint16_t>()
    { return "ushort"; }

template<> inline std::string get_nc_type<int32_t>()
    { return "int"; }

template<> inline std::string get_nc_type<uint32_t>()
    { return "uint"; }

template<> inline std::string get_nc_type<int64_t>()
    { return "int64"; }

template<> inline std::string get_nc_type<uint64_t>()
    { return "uint64"; }

template<> inline std::string get_nc_type<float>()
    { return "float"; }

template<> inline std::string get_nc_type<double>()
    { return "double"; }

template<> inline std::string get_nc_type<std::string>()
    { return "string"; }
// ---------------------------------------------------
/** Converts a string to a NetCDF type */
inline netCDF::NcType nc_type(netCDF::NcVar ncvar, std::string sntype)
    { return ncvar.getParentGroup().getType(sntype, netCDF::NcGroup::ParentsAndCurrent); }
// ---------------------------------------------------

struct TaggedThunk {
    std::function<void ()> fn;
    std::string tag;
    TaggedThunk(std::function<void ()> const &_fn, std::string const &_tag) :
        fn(_fn), tag(_tag) {}
};

/** Used to keep track of future writes on NcDefine */
class NcIO {
    std::vector<TaggedThunk> _io;

    netCDF::NcFile _mync;  // NcFile lacks proper move constructor
    bool own_nc;

    static void default_configure_var(netCDF::NcVar ncvar);

public:
    TmpAlloc tmp;    // Data kept around for the write phase
    netCDF::NcGroup * const nc;

    // Extracts of the read/write/define mode
    char const rw;
    const bool define;

    // Used to configure a variable after it's been created
    std::function<void(netCDF::NcVar)> const configure_var;

    // mode can be (see https://docs.python.org/3/library/functions.html#open)
    // 'r' 	open for reading (default)
    // 'w' 	open for writing, truncating the file first
    // 'x' 	open for exclusive creation, failing if the file already exists
    // 'a' 	open for writing, appending to the end of the file if it exists
    NcIO(std::string const &filePath, char mode,
        std::function<void(netCDF::NcVar)> const &_configure_var =
            std::bind(NcIO::default_configure_var, std::placeholders::_1));

    /** Create a "dummy" NcIO from an already-opened NetCDF file */
    NcIO(netCDF::NcGroup *_nc, char _rw) : own_nc(false), nc(_nc), rw(_rw), define(rw=='w') {}

    ~NcIO() { close(); }

    /** Converts a string to a NetCDF type */
    inline netCDF::NcType nc_type(std::string sntype)
        { return nc->getType(sntype, netCDF::NcGroup::ParentsAndCurrent); }


    void add(std::string const &tag, std::function<void ()> const &fn);

    void operator+=(std::function<void ()> const &fn)
        { add("", fn); }

    void flush(bool debug);
    void operator()()
        { flush(false); }

    void close();
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
    std::vector<long> const &dim_lens);

/** No unlimited dimensions possible here */
extern std::vector<netCDF::NcDim> get_or_add_finite_dims(
    NcIO &ncio,
    std::vector<std::string> const &dim_names,
    std::vector<size_t> const &dim_lens);

// ---------------------------------------------------------
template<class TypeT, int RANK>
std::vector<netCDF::NcDim> get_or_add_dims(
    NcIO &ncio,
    blitz::Array<TypeT, RANK> const &val,
    std::vector<std::string> const &dim_names);

template<class TypeT, int RANK>
std::vector<netCDF::NcDim> get_or_add_dims(
    NcIO &ncio,
    blitz::Array<TypeT, RANK> const &val,
    std::vector<std::string> const &dim_names)
{
    std::vector<long> dim_sizes(RANK);
    for (int k=0; k<RANK; ++k) dim_sizes[k] = val.extent(k);
    return get_or_add_dims(ncio, dim_names, dim_sizes);
}

// ---------------------------------------------------------
/** Accumulator for dimension name and sizes.
Used to concatenate dimensions from different places. */
class NcDimSpec {
    std::vector<std::string> names;
    std::vector<long> extents;

public:
    NcDimSpec() {}
    NcDimSpec(std::vector<std::string> &&_names,
        std::vector<long> &&_extents)
    : names(_names), extents(_extents) {}


    void push_back(std::string const &name, long extent) {
        names.push_back(name);
        extents.push_back(extent);
    }

    /** Convert to an array dimensions, which can be passed to
        get_or_add_var() */
    std::vector<netCDF::NcDim> to_dims(NcIO &ncio)
        { return get_or_add_dims(ncio, names, extents); }

};


// ===========================================================
// Variable Wrangling
netCDF::NcVar get_or_add_var(
    NcIO &ncio,
    std::string const &vname,
    std::string const &snc_type,
    std::vector<netCDF::NcDim> const &dims);

/*
void netCDF::NcVar::getVar 	( 	const std::vector< size_t > &  	start,
		const std::vector< size_t > &  	count,
		const std::vector< ptrdiff_t > &  	stride,
		const std::vector< ptrdiff_t > &  	imap,
		long long *  	dataValues 
	) 		const

Reads a mapped array section of values from a netCDF variable.

The mapped array section is specified by giving a corner, a vector of
edge lengths, a stride vector, and an index mapping vector. The index
mapping vector is a vector of integers that specifies the mapping
between the dimensions of a netCDF variable and the in-memory
structure of the internal data array. No assumptions are made about
the ordering or length of the dimensions of the data array.

Parameters

start Vector specifying the index in the variable where the first of
    the data values will be read. The indices are relative to 0, so
    for example, the first data value of a variable would have index
    (0, 0, ... , 0). The length of start must be the same as the
    number of dimensions of the specified variable. The elements of
    start correspond, in order, to the variable's dimensions. Hence,
    if the variable is a record variable, the first index would
    correspond to the starting record number for reading the data
    values.

count Vector specifying the edge lengths along each dimension of the
    block of data values to be read. To read a single value, for
    example, specify count as (1, 1, ... , 1). The length of count is
    the number of dimensions of the specified variable. The elements
    of count correspond, in order, to the variable's
    dimensions. Hence, if the variable is a record variable, the first
    element of count corresponds to a count of the number of records
    to read. Note: setting any element of the count array to zero
    causes the function to exit without error, and without doing
    anything.  stride Vector specifying the interval between selected
    indices. The elements of the stride vector correspond, in order,
    to the variable's dimensions. A value of 1 accesses adjacent
    values of the netCDF variable in the corresponding dimension; a
    value of 2 accesses every other value of the netCDF variable in
    the corresponding dimension; and so on. A NULL stride argument is
    treated as (1, 1, ... , 1).

imap Vector of integers that specifies the mapping between the
    dimensions of a netCDF variable and the in-memory structure of the
    internal data array. imap[0] gives the distance between elements
    of the internal array corresponding to the most slowly varying
    dimension of the netCDF variable. imap[n-1] (where n is the rank
    of the netCDF variable) gives the distance between elements of the
    internal array corresponding to the most rapidly varying dimension
    of the netCDF variable. Intervening imap elements correspond to
    other dimensions of the netCDF variable in the obvious
    way. Distances between elements are specified in type-independent
    units of elements (the distance between internal elements that
    occupy adjacent memory locations is 1 and not the element's
    byte-length as in netCDF 2).

dataValues Pointer to the location into which the data value is
    read. If the type of data value differs from the netCDF variable
    type, type conversion will occur. (However, no type conversion is
    carried out for variables using the user-defined data types:
    nc_Vlen, nc_Opaque, nc_Compound and nc_Enum.)

*/
template<class TypeT>
void get_or_put_var(netCDF::NcVar &ncvar, char rw,
    std::vector<size_t> const &start,
    std::vector<size_t> const &count,
    std::vector<ptrdiff_t> const &stride,
    std::vector<ptrdiff_t> const &imap,
    TypeT *dataValues);

template<class TypeT>
void get_or_put_var(netCDF::NcVar &ncvar, char rw,
    std::vector<size_t> const &start,
    std::vector<size_t> const &count,
    std::vector<ptrdiff_t> const &stride,
    std::vector<ptrdiff_t> const &imap,
    TypeT *dataValues)
{
    switch(rw) {
        case 'r' :
            ncvar.getVar(start, count, stride, imap, dataValues);
        break;
        case 'w' :
            ncvar.putVar(start, count, stride, imap, dataValues);
        break;
    }
}
// ========================================================
// Attribute Wrangling

template<class NcVarT, class AttrT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name, std::string const &stype,
    AttrT *data, size_t len);

template<class NcVarT, class AttrT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name, std::string const &stype,
    AttrT *data, size_t len)
{
    switch(rw) {
        case 'w': {
            if (std::is_same<AttrT,bool>::value) {
                std::string sbool;
                sbool.reserve(len);
                for (size_t i=0; i<len; ++i) sbool.push_back(data[i] ? 't' : 'f');
                ncvar.putAtt(name, sbool);
            } else {
                ncvar.putAtt(name, nc_type(ncvar, stype), len, data);
            }
        } break;
        case 'r': {
            auto att(ncvar.getAtt(name));
            if (att.getAttLength() != len) {
                (*ibmisc_error)(-1,
                    "Trying to read attribute %s of length %ld into C++ "
                    "variable of length %ld",
                    name.c_str(), att.getAttLength(), len);
            }

            if (std::is_same<AttrT,bool>::value) {
                std::string sbool;
                att.getValues(sbool);
                for (size_t i=0; i<len; ++i) {
                    data[i] = (
                        (sbool[i] == 't') || (sbool[i] == 'T') || (sbool[i] == '0'));
                }
            } else {
                att.getValues(data);
            }
        } break;
    }
}



template<class NcVarT, class AttrT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name,
    AttrT *data, size_t len);

template<class NcVarT, class AttrT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name,
    AttrT *data, size_t len)
{
    get_or_put_att(ncvar, rw, name, get_nc_type<AttrT>(), data, len);
}
// ---------------------------------------
template<class NcVarT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    std::string const &name,
    std::string &data,
    bool required = true);

template<class NcVarT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    std::string const &name,
    std::string &data,
    bool required)
{
    switch(rw) {
        case 'w':
            ncvar.putAtt(name, data);
        break;
        case 'r':
            auto att(ncvar.getAtt(name));
            if (required || !att.isNull())
                att.getValues(data);
        break;
    }
}

// ---------------------------------------
template<class NcVarT>
inline void get_or_put_att(
    NcVarT &ncvar, char rw,
    std::string const &name,
    bool &data)
{
    int idata;
    if (data) idata = 1;
    else idata = 0;
    get_or_put_att(ncvar, rw, name, "int", &idata, 1);
    if (rw == 'r') data = idata;
}
// ---------------------------------------
template<class NcVarT, class AttrT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name, std::string const &sntype,
    std::vector<AttrT> &data);

template<class NcVarT, class AttrT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name, std::string const &sntype,
    std::vector<AttrT> &data)
{
    auto type(nc_type(ncvar, sntype));
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

/** Overloaded function just for std::vector<std::string>, due to
    non-standard way NetCDF4-C++ handles that data type. */
template<class NcVarT>
extern void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name,
    std::string const &sntype,    // Ignored here...
    std::vector<std::string> &data);

template<class NcVarT>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name,
    std::string const &sntype,    // Ignored here...
    std::vector<std::string> &data)
{
    auto string_t(nc_type(ncvar, "string"));
    switch(rw) {
        case 'w':
            ncvar.putAtt(name, string_t, data.size(), &data[0]);
        break;
        case 'r':
            // String arrays are inexplicably returned as arrays
            // of char * that must be manually freed.  Convert
            // to civilized std::vector<std::string>
            auto att(ncvar.getAtt(name));
            auto N(att.getAttLength());
            std::vector<char *> cstrs(N);
            att.getValues(&cstrs[0]);
            data.clear();
            data.reserve(cstrs.size());
            for (size_t i=0; i<cstrs.size(); ++i) {
                data.push_back(std::string(cstrs[i]));
                free(cstrs[i]);
            }
        break;
    }
}

// ---------------------------------------
template<class NcVarT, class AttrT, size_t LEN>
void get_or_put_att(
    NcVarT &ncvar, char rw,
    const std::string &name, std::string const &sntype,
    std::array<AttrT,LEN> &data)
{
    auto data_v(to_vector(data));
    return get_or_put_att(ncvar, rw, name, sntype, data_v);
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
// ---------------------------------------

// ---------------------------------------


// ====================== Error Checking =============================
// ---------------------------------------------------
/** Check that blitz::Array is unit strides, row major.
For now, that's the only kind of Blitz variable we know how to write. */
template<class TypeT, int RANK>
void _check_blitz_strides(blitz::Array<TypeT, RANK> const &val, std::string const &vname);

template<class TypeT, int RANK>
void _check_blitz_strides(blitz::Array<TypeT, RANK> const &val, std::string const &vname)
{
    size_t expected_stride = 1;
    for (int k=RANK-1; k >= 0; --k) {
        size_t actual_stride = val.stride(k);
        if (val.stride(k) != expected_stride)
            (*ibmisc_error)(-1, "blitz::Array dimension %s-%d has unexpected stride %ld (vs %ld), cannot read/write with NetCDF.", vname.c_str(), k, actual_stride, expected_stride);
        expected_stride *= val.extent(k);
    }
}
// ---------------------------------------------------
/** Check that the NetCDF variable has the correct rank */
void _check_nc_rank(netCDF::NcVar const &ncvar, int rank);
// ---------------------------------------------------
template<class TypeT, int RANK>
void _check_blitz_dims(
    netCDF::NcVar const &ncvar,
    blitz::Array<TypeT, RANK> const &val,
    char rw);

template<class TypeT, int RANK>
void _check_blitz_dims(
    netCDF::NcVar const &ncvar,
    blitz::Array<TypeT, RANK> const &val,
    char rw)
{
    _check_nc_rank(ncvar, RANK);

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
void _check_vector_dims(
    netCDF::NcVar const &ncvar,
    std::vector<TypeT> const &val,
    char rw);

template<class TypeT>
void _check_vector_dims(
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
// =================================================
/** Typical array metadata found in NetCDF, etc. */
template<int RANK>
struct ArrayMeta {
    std::string name;
    std::array<int, RANK> shape;
    std::array<std::string,RANK> sdims;
    std::vector<std::tuple<std::string, std::string>> attr;    // (attr-name, value)

    /** Users do not use directly; see def() */
    ArrayMeta(
        std::string const &_name,
        std::array<int, RANK> const &_shape,
        std::array<std::string,RANK> _sdims,
        std::vector<std::tuple<std::string, std::string>> &&_attr)
    : name(_name), shape(_shape), sdims(_sdims), attr(std::move(_attr)) {}

public:
    /** Sets the shape of a bundle variable, but does not allocate. */
    void set_shape(
        std::array<int, RANK> const &_shape,
        std::array<std::string,RANK> const &_sdims,
        bool check = true)
    {
        if (check && shape[0] >= 0) (*ibmisc_error)(-1,
            "ArrayBundle variable %s shape already set", name.c_str());
        shape = _shape;
        sdims = _sdims;
    }

    template<class TypeT>
    void set_shape(
        blitz::Array<TypeT, RANK> const &arr,
        std::array<std::string,RANK> const &_sdims,
        bool check = true)
    {
        if (check && shape[0] >= 0) (*ibmisc_error)(-1,
            "ArrayBundle variable %s shape already set", name.c_str());

        for (int i=0; i<RANK; ++i) shape[i] = arr.extent(i);
        sdims = _sdims;
    }

    void set_shape(netCDF::NcVar &ncvar, int check=true)
    {
        auto dims(ncvar.getDims());
        if (dims.size() != RANK) (*ibmisc_error)(-1,
            "Rank mismatch (%ld vs %d) setting shape from NetCDF variable %s",
            dims.size(), RANK, ncvar.getName().c_str());

        for (int i=0; i<RANK; ++i) {
            netCDF::NcDim ncdim = ncvar.getDim(i);
            shape[i] = ncdim.getSize();
            sdims[i] = ncdim.getName();
        }
    }

#if 0
    /** Sets shape preferentially from (a) NetCDF variable, (b) Blitz array */
    template<class TypeT>
    void set_shape(
        netCDF::NcVar &ncvar,
        blitz::Array<TypeT, RANK> const &arr,
        std::array<std::string,RANK> const &_sdims,
        bool check = true)
    {
        // NcVar ncvar = nc->getVar(vname);
        if (!ncvar.isNull()) {
            // Set from NetCDF
            auto dims(ncvar.getDims());
            if (dims.size() != RANK) (*ibmisc_error)(-1,
                "Rank mismatch (%ld vs %d) setting shape from NetCDF variable %s",
                dims.size(), RANK, ncvar.getName().c_str());

            for (int i=0; i<RANK; ++i) {
                netCDF::NcDim ncdim = ncvar.getDim(i);
                std::string const ncdimname = ncdim.getName();
                shape[i] = ncdim.getSize();

                if (_sdims[i] != "" && _sdims[i] != ncdimname) (*ibmisc_error)(-1,
                    "Dimension name mismatch[%d]: %s vs %s",
                    i, _sdims[i].c_str(), ncdimname.c_str());
                sdims[i] = ncdim.getName();
            }
        } else {
            // Set from Blitz Array
            set_shape(arr, _sdims, check);
        }

    }
#endif

};    // struct Meta


// =================================================
enum class DimOrderMatch {
    MEMORY,        // Keep ordering of elements in memory the same
    LEXICAL };        // Keep ordering of dimensions the same

class NamedDim {
    std::string name;    // "" if unknown
    int extent;        // <0 if unknown
    NamedDim() : name(""), extent(-1) {}
};

std::vector<NamedDim> named_dims(netCDF::NcVar &ncvar)
{
    std::vector<NamedDim> ret;
    if (ncar.isNull()) return ret;

    for (int in=0; i<ncvar.getDimCount(); ++in) {
        NcDim ncdim(ncvar.getDim(in));
        ret.push_back(NamedDim(ncdim.getName(), ncdim.getSize()));
    }
    return ret;
}


std::vector<NamedDim> named_dims(
    std::vector<netCDF::NetCDF> const &ncdims,
    std::vector<int> const &ordering)
{
    bool const nc_dims_in_nc_order = (ordering.size() == 0);
    std::vector<NamedDim> ret;

    for (int in=0; i<ncvar.getDimCount(); ++in) {
        int const ib = (nc_dims_in_nc_order ? in : ordering.size()-in-1);
        netCDF::NcDim &ncdim(ncdims[ib]);

        if (ncdim.isNull()) {
            ret.push_back(NamedDim());    // slot / placeholder
        } else {
            ret.push_back(NamedDim(ncdim.getName(), ncdim.getSize()));
        }
    }
    return ret;
}

std::vector<NamedDim> named_dims(
    std::vector<netCDF::NetCDF> const &ncdims,
    std::vector<int> const &ordering,
    netCDF::NcVar &ncvar)
{
    std::vector<NamedDim> ncdims;
    if (dims.size() == 0) {
        ncdims = named_dims(ncvar);
    } else {
        ncdims = named_dims(dims, ordering);
    }
    return ncdims;
}


template<class TypeT, int RANK>
std::vector<NamedDim> named_dims(blitz::Array<TypeT, RANK> const &var)
{
    std::vector<NamedDim> ret;
    ret.reserve(RANK);
    if (var.data()) {    // Array is allocated
        for (int i=0; i<RANK; ++i) ret.push_back(NamedDim("", ret.extent(i)));
    } else {
        for (int i=0; i<RANK; ++i) ret.push_back(NamedDim());
    }
    return ret;
}

// ==============================================================


namespace _ncio_blitz {

template<class TypeT, int RANK>
void nc_rw_blitz2(
    netCDF::NcGroup *nc,
    char rw,
    blitz::Array<TypeT, RANK> *val,
    std::string const &vname,
    std::vector<size_t> const &nc_start,        // Where to the NetCDF variable; could have more than RANK dimensions
    std::array<int,RANK> const &b2n)    // (Blitz dim i) corresponds to (NetCDF dim b2n[i])
;

template<class TypeT, int RANK>
void nc_rw_blitz2(
    netCDF::NcGroup *nc,
    char rw,
    blitz::Array<TypeT, RANK> *val,
    std::string const &vname,
    std::vector<size_t> const &nc_start,        // Where to the NetCDF variable; could have more than RANK dimensions
    std::array<int,RANK> const &b2n)    // (Blitz dim i) corresponds to (NetCDF dim b2n[i])
{
    netCDF::NcVar ncvar = nc->getVar(vname);
    int const nc_rank = ncvar.getDimCount();

    // Set up vectors to be used in raw NetCDF call
    std::vector<size_t> start;
    std::vector<size_t> count;
    std::vector<ptrdiff_t> stride;
    std::vector<ptrdiff_t> imap;

    for (int i=0; i<nc_rank; ++i) {
        start.push_back(nc_start[i]);
        count.push_back(1);
        stride.push_back(1);    // Stride in NetCDF variable; maybe let the user set this?
        imap.push_back(1);        // No stride for dimensions that don't correspond to Blitz var
    }

    for (int i=0; i<RANK; ++i) {
        imap[b2n[i]] = val->stride(i);
        count[b2n[i]] = val->extent(i);
        // Allow nc_dim[b2n[i]] > val->extent(i)
    }

    get_or_put_var(ncvar, rw, start, count, stride, imap, val->data());
}
// -------------------------------------------------------------------



/** Parameters passed to ncio_blitz() by helper functions */
struct Info {
    // Output of
    std::vector<NamedDim> blitz, netcdf;
    std::vector<int> b2n;  // (Blitz dim i) corresponds to (NetCDF dim b2n[i]).  If {}, then {0,1,2,...}; otherise, rank=RANK

    // Info about NetCDF variable
    std::vector<netCDF::NcDim> dims;    // Dimensions of ON-DISK NetCDF Variable in NetCDF order
    std::vector<size_t> nc_start;        // Where to the NetCDF variable; same rank as dims.   If {}, then all 0's.  <0 means "allocate a Blitz var for this dim."


    /** Reconciles dimensions between Blitz and NetCDF.  This subroutine
        fills in the blanks, it does not check afterwards.  That is for elsewhere. */
    Info::Info(
        std::vector<NamedDim> &&blitz,
        std::vector<int> const &blitz_ordering,    // Ordering of dimensions in increasing stride
        DimOrderMatch match,
        std::vector<NamedDim> &&netcdf,
        // bool ncdims_in_nc_order,    == true ALWAYS
        std::vector<int> &&b2n,
        std::vector<size_t> &&nc_start);

    

};

// ------------------------------------------------------------------------------
template<class TypeT, int RANK>
struct Helper {
public:
    std::vector<netCDF::NcDim> dims;

    Helper(std::vector<netCDF::NcDim> const &_dims) : dims(_dims) {}
    virtual ~Helper() {}

    virtual std::vector<int> ordering() = 0;

    virtual Info<RANK> _help(
        NcIO &ncio, netCDF::NcVar &ncvar, blitz::Array<TypeT,RANK> &arr) = 0;

    Info<RANK> operator()(
        NcIO &ncio, std::string const &vname, blitz::Array<TypeT,RANK> &arr)
    {
        // nc->getVar() returns a null NcVar if no object of that name is found.
        // That is what we want here.
        netCDF::NcVar ncvar = ncio.nc->getVar(vname);
        return _help(ncio, ncvar, arr);
    }

};

#define NCIO_BLITZ_PARAMS \
    NcIO &ncio, \
    blitz::Array<TypeT, RANK> &arr, \    // Allocated
    std::string const &vname, \
    std::string const &snc_type
#define NCIO_BLITZ_ARGS \
    ncio, arr, vname, snc_type

/** Maps a blitz::Array directly to/from a NetCDF array
@param val The array to read/write
@param vname Name of NetCDF variable to write
@param snc_type Data type to write in NetCDF variable.  See get_nc_type().
@param dims Dimensions to define NetCDF variable.  See get_or_add_dims(). */
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    Info<RANK> const &info);

/** Maps a blitz::Array directly to/from a NetCDF array */
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    Info<RANK> const &info)
{
    netCDF::NcVar ncvar(ncio.nc->getVar(vname));

    if (info.blitz_rank() != RANK) (*ibmisc_error)(-1,
        "Rank of Info=%d must match Blitz rank=%d",
        info.blitz_rank(), RANK);

    // Get NetCDF variable
    int const nc_rank = info.dims.size();
    if (ncio.rw == 'w') {
        // Replace possibly null NcVar with one we've checked out and maybe created.
        ncvar = get_or_add_var(ncio, vname,
            snc_type == "" ? get_nc_type<TypeT>() : snc_type,
            info.dims);
    }
printf("BB4\n");

    // Check RANK vs. info.dims.size()
    if (RANK > info.dims.size()) (*ibmisc_error)(-1,
        "Formal NetCDF RANK=%ld must be AT LEAST Blitz RANK=%d",
        info.dims.size(), RANK);
printf("BB5\n");

    // Check info.dims vs. dimensions on disk
    // NetCDF4-C++ library does not bounds check (as of 2016-01-15)
    if (ncvar.getDimCount() != info.dims.size()) (*ibmisc_error)(-1,
        "Existing NetCDF RANK=%d must equal formal rank=%ld",
        ncvar.getDimCount(), info.dims.size());
printf("BB6\n");

    for (int k=0; k<info.dims.size(); ++k) {
        // The dimensions must match, not just their extents
        if (info.dims[k] != ncvar.getDim(k)) (*ibmisc_error)(-1,
            "User-supplied dimension[%d] does not match netCDF dimension=%s",
            k, ncvar.getDim(k).getName().c_str());
    }
printf("BB7\n");

    // Check: in-memory variable is compatible with NetCDF variable
    for (int kb=0; kb<RANK; ++kb) {
        int const kn = info.b2n[kb];
        if (arr.extent(kb) > ncvar.getDim(kn).getSize() - info.nc_start[kn]) (*ibmisc_error)(-1,
            "C++ variable extent[%d]=%d does not match NetCDF variable extent[%d]=%d",
            kb, arr.extent(kb), kn, ncvar.getDim(kn).getSize(), info.nc_start[kn]);
    }

    // const_cast allows us to re-use nc_rw_blitz for read and write
    ncio += std::bind(&nc_rw_blitz2<TypeT, RANK>,
        ncio.nc, ncio.rw, &arr, vname, info.nc_start, info.b2n);

    return ncvar;
}
// -----------------------------------------------------------------------------

// ------------------------------------------------------------------------------
/** read/write entire NetCDF variable into an existing Blitz array
@param blitz_order
    true: dim_names/dim_lens correspond to order of dimensions in arr.
          Order of dimensions in NetCDF var. will be a row-major reordering
          of these dimensions.  This results in the same ordering of data
          on disk as in memory.
    false: dim_names/dim_lens correspond to order of dimensions in NetCDF variable.
          This could result in different ordering of data on disk vs. memory.

@param equal_dim_order
    true: Dimensions of Blitz++ and NetCDF variable are in the same
          lexical order, even if strides differ.

    false: Blitz++ and NetCDF variable are laid out the same in
          memory, even if their lexical dimension ordering differs due
          to, eg. column-major vs. row-major.
*/
template<class TypeT, int RANK>
struct Helper_whole1 : public Helper<TypeT,RANK>
{
    bool equal_dim_order;
    bool dims_in_nc_order;

    Helper_whole1(
        std::vector<netCDF::NcDim> const &_dims,
        DimOrderMatch _match,
        bool _dims_in_nc_order)
    : Helper<TypeT,RANK>(_dims), match(_match),
    dims_in_nc_order(_dims_in_nc_order) {}

    Info _help(
        // Leave these unbound
        NcIO &ncio,
        netCDF::NcVar &ncvar,
        blitz::Array<TypeT, RANK> &arr)
    {
printf("BEGIN _whole1\n");
        auto &dims(Helper<TypeT,RANK>::dims);    // Template problem accessing superclass member
        int const nc_rank = RANK;

        if (!arr.data()) (*ibmisc_error)(-1,
            "Blitz array must be pre-allocated for Helper_whole1; try using Helper_alloc");

        auto ordering(to_vector(arr.ordering()));
        auto ncdims(named_dims(dims, dims_in_nc_order ? {} : ordering, ncvar));
        Info ret(
            named_dims(arr), ordering,
            match,
            std::move(ncdims), ncvar),
            {}, {});    // b2n, nc_start

        return ret;
    }
};
// -------------------------------------------------------------------------------
template<class TypeT, int RANK>
struct Helper_alloc1 : public Helper<TypeT,RANK>
{
    bool dims_in_nc_order;
    blitz::GeneralArrayStorage<RANK> storage;

    Helper_alloc1(
        std::vector<netCDF::NcDim> const &_dims,
        bool _dims_in_nc_order,
        blitz::GeneralArrayStorage<RANK> const &_storage)
    : Helper<TypeT,RANK>(std::move(_dims)),
    dims_in_nc_order(_dims_in_nc_order),
    storage(_storage) {}

    Info<RANK> _help(
        // Leave these unbound
        NcIO &ncio,
        netCDF::NcVar &ncvar,
        blitz::Array<TypeT, RANK> &arr)
    {
        auto &dims(Helper<TypeT,RANK>::dims);    // Template problem accessing superclass member

        // Array in memory must (will be allocated to) have same
        // memory layout as on disk.
        DimOrderMatch const match = DimOrderMatch::MEMORY;


        auto ordering(to_vector(storage.ordering()));
        auto ncdims(named_dims(dims, dims_in_nc_order ? {} : ordering, ncvar));
        Info ret(
            named_dims(arr), ordering,
            match,
            std::move(ncdims), ncvar),
            {}, {});    // b2n, nc_start

        // Allocate the Blitz array
        // This only makes sense for reading
        if (ncio.rw == 'r') {

            blitz::TinyVector<int,RANK> extent;
            for (int ib=0; ib<RANK; ++ib) extent[ib] = ret.blitz[ib].extent;
            arr.reference(blitz::Array<TypeT,RANK>(extent, storage));
        }

        return ret;
    }
};
// -------------------------------------------------------------------------
template<class TypeT, int RANK>
struct Helper_partial : public Helper<TypeT,RANK>
{
    std::vector<int> nc_start;    // Where to start each dimension in NetCDF
    std::array<int,RANK> const &b2n;    // Where to slot each Blitz++ dimension

    Helper_partial(
        std::vector<netCDF::NcDim> const &_dims,
        std::vector<int> _nc_start,    // Where to start each dimension in NetCDF
        std::array<int,RANK> const &_b2n)    // Where to slot each Blitz++ dimension
    : Helper<TypeT,RANK>(std::move(_dims)),
    nc_start(_nc_start),
    b2n(_b2n) {}


    Info<RANK> _help(
        // Leave these unbound
        NcIO &ncio,
        netCDF::NcVar &ncvar,
        blitz::Array<TypeT, RANK> &arr)
    {
        auto &dims(Helper<TypeT,RANK>::dims);    // Template problem accessing superclass member

        // By necessity, we must list dimensions in NetCDF order, so we know
        // where Blitz++ dimensions slot in
        bool const dims_in_nc_order=true;

        auto ordering(to_vector(storage.ordering()));
        auto ncdims(named_dims(dims, dims_in_nc_order ? {} : ordering, ncvar));

        Info ret(
            named_dims(arr), arr.ordering(),
            DimMatch::MEM,    // Not used
            std::move(ncdims),
            b2n, nc_start);

            named_dims(dims, dims_in_nc_order ? {} : ordering, n


        return ret;
    }
}    ;    
}    // namespace _ncio_blitz
// ==================================================================
template<class TypeT, int RANK>
inline netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &dims,
    bool equal_dim_order=false,
    bool dims_in_nc_order=true)
{
    using namespace std::placeholders;

    _ncio_blitz::Helper_whole1<TypeT,RANK> helper(
            dims, equal_dim_order, dims_in_nc_order);
    return _ncio_blitz::ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS,
        helper(ncio, vname, arr));
}

template<class TypeT, int RANK>
inline netCDF::NcVar ncio_blitz_alloc(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &dims = {},
    bool dims_in_nc_order = true,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>())
{
    using namespace std::placeholders;
    _ncio_blitz::Helper_alloc1<TypeT,RANK> helper(
        dims, dims_in_nc_order, storage);
    return _ncio_blitz::ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS,
        helper(ncio, vname, arr));
}

template<class TypeT, int RANK>
inline netCDF::NcVar ncio_blitz_partial(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &dims,
    std::vector<int> const &nc_start,    // Where to start each dimension in NetCDF
    std::array<int,RANK> const &b2n)    // Where to slot each Blitz++ dimension
{
    using namespace std::placeholders;
    _ncio_blitz::Helper_partial<TypeT,RANK> helper(
        dims, nc_start, b2n);
    return _ncio_blitz::ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS,
        helper(ncio, vname, arr));
}
#undef NCIO_BLITZ_PARAMS


/** Convenience method: reads data from NetCDF, returns as a newly allocate blitz::Array */
template<class TypeT, int RANK>
inline blitz::Array<TypeT,RANK> nc_read_blitz(
    netCDF::NcGroup *nc,
    std::string const &vname,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>())
{
    NcIO ncio(nc, 'r');    // Dummy
    blitz::Array<TypeT,RANK> val;
    ncio_blitz_alloc(ncio, val, vname, get_nc_type<TypeT>(), {}, true, storage);
    return val;
}
// =================================================
// Specializations for std::vector instead of blitz::Array

/** Like get_or_add_dims() above, but specialized for std::vector */
template<class TypeT>
std::vector<netCDF::NcDim> get_or_add_dims(
    NcIO &ncio,
    std::vector<TypeT> &val,
    std::array<std::string, 1> const &dim_names)
{
    std::array<size_t, 1> dim_sizes;
    dim_sizes[0] = val.size();
    return get_or_add_dims(ncio, val, dim_names, dim_sizes);
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
    _check_nc_rank(ncvar, 1);

    size_t ncsize = ncvar.getDim(0).getSize();

    if (rw == 'r') {
        if (alloc) val->resize(ncsize);
        // Vector dimensions checked below
    }

    _check_vector_dims(ncvar, *val, rw);

    std::vector<size_t> startp = {0};
    std::vector<size_t> countp = {ncsize};
    switch(rw) {
        case 'r' :
            ncvar.getVar(startp, countp, &(*val)[0]);
        break;
        case 'w' :
            ncvar.putVar(startp, countp, &(*val)[0]);
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
    bool alloc,         // Should we allocate val?
    std::string const &vname,
    std::string const &snc_type,
    std::vector<netCDF::NcDim> const &dims)
{
    get_or_add_var(ncio, vname, snc_type, dims);

    ncio += std::bind(&nc_rw_vector<TypeT>, ncio.nc, ncio.rw, &val, alloc, vname);
}
// ----------------------------------------------------

// ----------------------------------------------------

/** Do linewrap for strings that are intended to be used as comment attributes in NetCDF files.
       see: http://www.cplusplus.com/forum/beginner/19034/
*/
extern std::string ncwrap( std::string const &str, size_t width = 55 );


}   // Namespace
#endif  // Guard
