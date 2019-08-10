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
#include <ncException.h>
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
netCDF::NcType nc_type(netCDF::NcVar const &ncvar, std::string sntype);
netCDF::NcType nc_type(netCDF::NcGroup const &ncgroup, std::string sntype);
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

    std::unique_ptr<netCDF::NcFile> _mync;  // NcFile lacks proper move constructor

public:
    static void default_configure_var(netCDF::NcVar ncvar);
    static void no_compress(netCDF::NcVar ncvar);

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
    NcIO(std::string const &filePath, char mode = 'r',
        std::string const &format = "nc4",
        std::function<void(netCDF::NcVar)> const &_configure_var =
            std::bind(NcIO::default_configure_var, std::placeholders::_1));

    /** Create a "dummy" NcIO from an already-opened NetCDF file */
    NcIO(netCDF::NcGroup *_nc, char _rw) : nc(_nc), rw(_rw), define(rw=='w') {}

    ~NcIO() { close(); }

    /** Converts a string to a NetCDF type */
    inline netCDF::NcType nc_type(std::string sntype)
        { return nc->getType(sntype, netCDF::NcGroup::ParentsAndCurrent); }


    void add(std::string const &tag, std::function<void ()> const &fn);

    void operator+=(std::function<void ()> const &fn)
        { add("", fn); }

    void flush(bool debug=false);
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
@param dim_size Size to create or check for.  NOTES:
       1. If <0, assume user is requesting an unlimited dimension.
       2. Otherwise, this is checked when writing, ignored when reading.
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

/** Runs get_or_add_dim() on a number of dimensions
@param ncio The NetCDF file to read/write
@param dim_names Names of dimension to create.
@param dim_sizes Sizes to create or check for.  NOTES:
       1. If <0, assume user is requesting an unlimited dimension.
       2. Otherwise, this is checked when writing, ignored when reading.
       3. Must be same length as dim_names.
@return The created/retrieved dimensions.
*/
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

/* From NetCDF-cxx4 Manual:

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
    anything.

stridep A vector of ptrdiff_t integers that specifies the sampling
    interval along each dimension of the netCDF variable. The elements
    of the stride vector correspond, in order, to the netCDF
    variable's dimensions (stride[0] gives the sampling interval along
    the most slowly varying dimension of the netCDF
    variable). Sampling intervals are specified in type-independent
    units of elements (a value of 1 selects consecutive elements of
    the netCDF variable along the corresponding dimension, a value of
    2 selects every other element, etc.). A NULL stride argument is
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
    // NetCDF library is really slow if stride or count is used.
    // See if we can do without...
    ptrdiff_t cur_imap = 1;
    for (int i=imap.size()-1; i >= 0; --i) {
        if (imap[i] != cur_imap) {
            // Non-contiguous array; must do with full imap
            switch(rw) {
                case 'r' :
                    ncvar.getVar(start, count, stride, imap, dataValues);
                break;
                case 'w' :
                    ncvar.putVar(start, count, stride, imap, dataValues);
                break;
            }
            return;
        }
        cur_imap *= count[i];
    }

    // Contiguous array; can omit stride and imap
    switch(rw) {
        case 'r' :
            ncvar.getVar(start, count, dataValues);
        break;
        case 'w' :
            ncvar.putVar(start, count, dataValues);
        break;
    }
}
// ========================================================

// Attribute Wrangling

extern netCDF::NcVarAtt get_att(netCDF::NcVar &ncvar, std::string const &name);

extern netCDF::NcGroupAtt get_att(netCDF::NcGroup &ncvar, std::string const &name);

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
            auto att(get_att(ncvar, name));
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
            auto att(get_att(ncvar, name));

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
            auto att(get_att(ncvar, name));
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
            auto att(get_att(ncvar, name));
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
    if (rw == 'w') {
        auto data_v(to_vector(data));
        get_or_put_att(ncvar, rw, name, sntype, data_v);
    } else {
        std::vector<AttrT> data_v;
        get_or_put_att(ncvar, rw, name, sntype, data_v);
        if (data_v.size() != LEN) (*ibmisc_error)(-1,
            "Attribute %s has wrong length: %d instead of %d",
            name.c_str(), (int)data_v.size(), (int)LEN);

        for (int i=0; i<LEN; ++i) data[i] = data_v[i];
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
            auto att(get_att(ncvar, name));
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
            "Size (%ld; var=%p base=%p) of std::vector must match %s:%s (%ld) in NetCDF",
            val.size(), &val, &val[0], ncvar.getName().c_str(),
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
    std::vector<std::pair<std::string, std::string>> attr;    // (attr-name, value)

    /** Users do not use directly; see def() */
    ArrayMeta(
        std::string const &_name,
        std::array<int, RANK> const &_shape,
        std::array<std::string,RANK> _sdims,
        std::vector<std::pair<std::string, std::string>> &&_attr)
    : name(_name), shape(_shape), sdims(_sdims), attr(std::move(_attr)) {}

public:
    /** Sets the shape of a bundle variable, but does not allocate. */
    void set_shape(
        std::array<int, RANK> const &_shape,
        std::array<std::string,RANK> const &_sdims,
        bool check = true)
    {
        if (check && shape[0] >= 0) (*ibmisc_error)(-1,
            "Bundle variable %s shape already set", name.c_str());
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
            "Bundle variable %s shape already set", name.c_str());

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
    std::map<std::string, std::string> make_attr_map() const
        { return std::map<std::string, std::string>(attr.begin(), attr.end()); }
};    // struct Meta


// =================================================
/** Argument explaining how ncio_blitz_*() should match dimensions between
the Blitz++ and NetCDF arrays. */
enum class DimOrderMatch {
    MEMORY,        // Keep ordering of elements in memory the same
    LEXICAL };        // Keep ordering of dimensions the same

/** The name and extent for a dimension.  This information may be abstracted
from any of:
  a) Blitz++ array (no names available)
  b) Blitz++ array plus dimension string labels
  c) NetCDF array
  d) User spec. for NetCDF array

Once dimension information is obtained in a NamedDim, it may
be manipulated and reconciled the same way, no matter where it
came from.
*/
struct NamedDim {
    std::string name;    // "" if unknown
    int extent;        // <0 if unknown
    NamedDim() : name(""), extent(-1) {}
    NamedDim(std::string const &_name, int _extent)
        : name(_name), extent(_extent) {}
};

/** Extract NamedDims from an existing NetCDF variable.
  a) Dimensions will be returned in NetCDF order.
  b) If the NetCDF variable does not exist (eg ncvar.isNull()),
     then a zero-length vector will be returned. */
extern std::vector<NamedDim> named_dims(netCDF::NcVar &ncvar);

/** Extract NamedDims from a set of dimensions provided by the user.
  a) Dimensions will be returned in NetCDF order.
  
@param ordering Describes the ordering of the dimensions provided by
  the user, with the dimension of smallest stride first.
  If {}, then NetCDF order (row major) is assumed. */
extern std::vector<NamedDim> named_dims(
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<int> const &ordering);

/** Extracts NamedDims from, in descending priority:
  a) User-provided dimensions (ncdims, ordering) if they exist
  b) NetCDF variable, if it exists.
*/
extern std::vector<NamedDim> named_dims(
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<int> const &ordering,
    netCDF::NcVar &ncvar);

/** Extracts NamedDims from, in descending priority:
  a) User-provided dimensions (ncdims, dims_in_nc_order, blitz_ordering) if they exist
  b) NetCDF variable, if it exists.
@param dims_in_nc_order true if the dims variable are in NetCDF order
   (false if they are in the order specified by ordering)
*/
extern std::vector<NamedDim> named_dims(
    std::vector<netCDF::NcDim> const &dims,
    bool dims_in_nc_order,
    std::vector<int> const &blitz_ordering,
    netCDF::NcVar &ncvar);

/** Extracts NamedDims from a Blitz variable.
If dimension names are not provided, they will be given the name ""
@param var The Blitz array providing the extents.
@sdims (optional) Names of the dimensions. */
template<class TypeT, int RANK>
std::vector<NamedDim> named_dims(blitz::Array<TypeT, RANK> const &var,
    std::vector<std::string> const &sdims={});
    
template<class TypeT, int RANK>
std::vector<NamedDim> named_dims(blitz::Array<TypeT, RANK> const &var,
    std::vector<std::string> const &sdims)
{
    std::vector<NamedDim> ret;
    ret.reserve(RANK);
    std::string const blank("");
    for (int i=0; i<RANK; ++i) {
        std::string const &sdim(sdims.size() == 0 ? blank : sdims[i]);
        ret.push_back(var.data() ?    // Is array allocated?
            NamedDim(sdim, var.extent(i))
            : NamedDim(sdim, -1));
    }
    return ret;
}

// ==============================================================

// Namespace containing ncio_blitz() helper stuff
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
        int i_n = b2n[i];
        if (i_n < 0 || i_n >= imap.size()) (*ibmisc_error)(-1,
            "b2n[%d]=%d is out of range (0--%ld)",
            i, b2n[i], imap.size());

        imap[b2n[i]] = val->stride(i);
        count[b2n[i]] = val->extent(i);
        // Allow nc_dim[b2n[i]] > val->extent(i)
    }


    get_or_put_var(ncvar, rw, start, count, stride, imap, val->data());
}
// -------------------------------------------------------------------



/** Result of reconciling NamedDims from Blitz and NetCDF variables.
All items will exist, even if the underlying variable does not exist.
(In this case, information is taken from the other variable) */
struct Info {
    /** Dimensions for Blitz array, in its native order (blitz::Array.ordering()) */
    std::vector<NamedDim> blitz;

    /** Dimensions for NetCDF array, in NetCDF (row major).
    NOTE: NetCDF rank may be larger than Blitz rank. */
    std::vector<NamedDim> netcdf;

    /** Converts between Blitz++ and NetCDF dimensions.
          in = b2n(ib)
    where:
        in = NetCDF dimension index (0...NetCDF rank)
        ib = Blitz dimension index (0...Blitz rank)
    */
    std::vector<int> b2n;  // (Blitz dim i) corresponds to (NetCDF dim b2n[i]).  If {}, then {0,1,2,...}; otherise, rank=RANK

    // Info about NetCDF variable
    /** NetCDF dimensions of the on-disk NetCDF variable, in NetCDF order. */
    std::vector<netCDF::NcDim> dims;
    /** Where to start rewding/writing the NetCDF variable; NetCDF rank. */
    std::vector<size_t> nc_start;


    /** Reconciles dimensions between Blitz and NetCDF.
    @param blitz NamedDims vector extracted from Blitz array,
        Set to {} if Blitz variable is to be allocated based on NetCDF variable.
        (see named_dims() functions above)
    @param blitz_ordering Order of dimensions in blitz parameter.
         Obtained either from blitz::Array::ordering() or blitz::StorageOrder::ordering()
    @param match Details how to match Blitz and NetCDF dimensions, if b2n is not set.
         (Only used where it makes sense.  If possible, Blitz and NetCDF dimensions are
         matched by string matching).
    @param netcdf NamedDims vector extracted from NetCDF variable, in NetCDF order.
        NOTES:
          a) Set to {} if NetCDF variable is to be allocated based on Blitz variable.
          b) see named_dims() functions above
          c) Dimensions set to NamedDim() are placeholders; and can be "filled in" by Blitz dimensions.
    @param b2n Mapping from Blitz to NetCDF dimension indices.
        If not user-specified, leave blank; will be determined automatically by matching
        dimensions.
    @param nc_start Position in each NetCDF dimension to start reading/writing.
    NOTES:
       1. Does not check whether the NamedDim arrays correspond to what is
          actually on disk or in the Blitz array.  Those checks are done later.
    */
    Info(
        std::vector<NamedDim> &&blitz,
        std::vector<int> const &blitz_ordering,    // Ordering of dimensions in increasing stride
        DimOrderMatch match,
        std::vector<NamedDim> &&netcdf,
        std::vector<int> &&b2n,
        std::vector<size_t> &&nc_start);

    int blitz_rank() const { return blitz.size(); }
};

// ------------------------------------------------------------------------------

#define NCIO_BLITZ_PARAMS \
    NcIO &ncio, \
    blitz::Array<TypeT, RANK> &arr, \
    std::string const &vname, \
    std::string const &snc_type
#define NCIO_BLITZ_ARGS \
    ncio, arr, vname, snc_type

/** Reads/writes a blitz::Array directly to/from a NetCDF array
@param ncio The open NetCDF file to read/write.
@param arr The Blitz array to read/write
@param vname Name of NetCDF variable to write
@param snc_type Data type to write in NetCDF variable.  See get_nc_type().
@param info More details on how to map Blitz and NetCDF variables. */
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(NCIO_BLITZ_PARAMS,
    Info &info);

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    Info &info)
{
    // Set up info.dims
    info.dims.clear();
    for (int in=0; in<info.netcdf.size(); ++in) {
        auto const &nd(info.netcdf[in]);
        info.dims.push_back(
            nd.extent < 0 ?
                // We don't know what this dim should be; hopefully it exists on disk
                ncio.nc->getDim(nd.name)
                // We do know what the dim should be; create or get/verify it
                : get_or_add_dim(ncio, nd.name, nd.extent));
    }

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

    // Check RANK vs. info.dims.size()
    if (RANK > info.dims.size()) (*ibmisc_error)(-1,
        "Formal NetCDF RANK=%ld must be AT LEAST Blitz RANK=%d",
        info.dims.size(), RANK);

    // Check info.dims vs. dimensions on disk
    // NetCDF4-C++ library does not bounds check (as of 2016-01-15)
    if (ncvar.getDimCount() != info.dims.size()) (*ibmisc_error)(-1,
        "Existing NetCDF RANK=%d must equal formal rank=%ld",
        ncvar.getDimCount(), info.dims.size());

    for (int k=0; k<info.dims.size(); ++k) {
        // The dimensions must match, not just their extents
        if (info.dims[k] != ncvar.getDim(k)) (*ibmisc_error)(-1,
            "User-supplied dimension[%d]=%s does not match netCDF dimension=%s",
            k, info.dims[k].getName().c_str(), ncvar.getDim(k).getName().c_str());
    }

    // Check: in-memory variable is compatible with NetCDF variable
    for (int kb=0; kb<RANK; ++kb) {
        int const kn = info.b2n[kb];
        if (arr.extent(kb) > ncvar.getDim(kn).getSize() - info.nc_start[kn]) (*ibmisc_error)(-1,
            "C++ variable extent[%d]=%d does not match NetCDF variable extent[%d]=%d",
            kb, arr.extent(kb), kn, ncvar.getDim(kn).getSize(), info.nc_start[kn]);
    }

    // const_cast allows us to re-use nc_rw_blitz for read and write
    ncio += std::bind(&nc_rw_blitz2<TypeT, RANK>,
        ncio.nc, ncio.rw, &arr, vname, info.nc_start, to_array<int,int,RANK>(info.b2n));

    return ncvar;
}
// -----------------------------------------------------------------------------
}    // namespace _ncio_blitz
// ==================================================================


/** Reads/writes a blitz::Array directly to/from a whole NetCDF array
@param ncio The open NetCDF file to read/write.
@param arr The Blitz array to read/write
@param vname Name of NetCDF variable to write
@param snc_type Data type to write in NetCDF variable.  See get_nc_type().
@param info More details on how to map Blitz and NetCDF variables.
@param (OPTIONAL) ncdims Dimensions to use for the NetCDF variable.
    If these are to be deduced, use {}.
@param match
    if MEMORY: Use same order of elements on disk as in memory, even if dimensions
         come out in different order.
    if LEXICAL: Use same order of dimensions in Blitz and NetCDF variables,
         even if this causes the order of elements to change.
@param ncdims_in_nc_order True if ncdims argument lists dimensions in the order on disk.
    Else dimensions are listed in same order as in Blitz array.
@param arr_sdims (OPTIONAL) Names of each Blitz dimension.
    If name metadata is not available, set to {}.
*/
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    DimOrderMatch match=DimOrderMatch::MEMORY,
    bool ncdims_in_nc_order=true,
    std::vector<std::string> const &arr_sdims = {});

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    DimOrderMatch match,
    bool ncdims_in_nc_order,
    std::vector<std::string> const &arr_sdims)
{
    netCDF::NcVar ncvar = ncio.nc->getVar(vname);

    if (ncio.rw == 'r' && ncvar.isNull()) (*ibmisc_error)(-1,
        "Trying to read from non-existant NetCDF variable %s", vname.c_str());

    if (!arr.data()) (*ibmisc_error)(-1,
        "Blitz array must be pre-allocated for ncio_blitz(); try using ncio_blitz_alloc() instead");

    auto ordering(tiny_to_vector(arr.ordering()));
    auto ncdims_nd(named_dims(ncdims, ncdims_in_nc_order, ordering, ncvar));    // Output must be in nc order
    _ncio_blitz::Info info(
        named_dims(arr, arr_sdims), ordering,
        match,
        std::move(ncdims_nd),
        {}, {});    // b2n, nc_start

    return _ncio_blitz::ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS, info);
}

/** Reads/writes a blitz::Array directly to/from a whole NetCDF array.
Allocates the blitz::Array if reading.
@param ncio The open NetCDF file to read/write.
@param arr The Blitz array to read/write
@param vname Name of NetCDF variable to write
@param snc_type Data type to write in NetCDF variable.  See get_nc_type().
@param info More details on how to map Blitz and NetCDF variables.
@param (OPTIONAL) ncdims Dimensions to use for the NetCDF variable.
    If these are to be deduced, use {}.
@param match
    if MEMORY: Use same order of elements on disk as in memory, even if dimensions
         come out in different order.
    if LEXICAL: Use same order of dimensions in Blitz and NetCDF variables,
         even if this causes the order of elements to change.
@param ncdims_in_nc_order True if ncdims argument lists dimensions in the order on disk.
    Else dimensions are listed in same order as in Blitz array.
@param storage Storage scheme used to allocate Blitz array.
@param arr_sdims (OPTIONAL) Names of each Blitz dimension.
    If name metadata is not available, set to {}.
*/
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_alloc(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims = {},
    DimOrderMatch match=DimOrderMatch::MEMORY,
    bool ncdims_in_nc_order = true,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>(),
    std::vector<std::string> const &arr_sdims = {});

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_alloc(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    DimOrderMatch match,
    bool ncdims_in_nc_order,
    blitz::GeneralArrayStorage<RANK> const &storage,
    std::vector<std::string> const &arr_sdims)
{
    netCDF::NcVar ncvar = ncio.nc->getVar(vname);
//    // Array in memory must (will be allocated to) have same
//    // memory layout as on disk.
//    DimOrderMatch const match = DimOrderMatch::MEMORY;

    auto ordering(tiny_to_vector(arr.data() ? arr.ordering() : storage.ordering()));
    auto ncdims_nd(named_dims(ncdims, ncdims_in_nc_order, ordering, ncvar));

    _ncio_blitz::Info info(
        named_dims(arr, arr_sdims), ordering,
        match,
        std::move(ncdims_nd),
        {}, {});    // b2n, nc_start

    // Allocate the Blitz array
    // This only makes sense for reading
    if (ncio.rw == 'r') {

        blitz::TinyVector<int,RANK> extent;
        for (int ib=0; ib<RANK; ++ib) extent[ib] = info.blitz[ib].extent;
        arr.reference(blitz::Array<TypeT,RANK>(extent, storage));
    }


    return _ncio_blitz::ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS, info);
}

/** Reads/writes a blitz::Array directly to/from part of a NetCDF array.
@param ncio The open NetCDF file to read/write.
@param arr The Blitz array to read/write
@param vname Name of NetCDF variable to write
@param snc_type Data type to write in NetCDF variable.  See get_nc_type().
@param ncdims Dimensions to use for the NetCDF variable.
    If placeholders are to be used to slot Blitz dimensions into NetCDF variable, use NamedDim().
    If these are to be read from existing NetCDF variable, use {}.
@param nc_start Position in each NetCDF dimension to start reading/writing.
@param b2n (OPTIONAL) Mapping from Blitz to NetCDF dimension indices.
@param arr_sdims (OPTIONAL) Name of Blitz dimensions.
    If names are available for Blitz and NetCDF dimensions, they will be mapped
    by name, rather than relying on user-provided b2n. */
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_partial(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<size_t> const &nc_start,    // Where to start each dimension in NetCDF
    std::vector<int> const &b2n,    // Where to slot each Blitz++ dimension
    std::vector<std::string> const &arr_sdims = {});

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_partial(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<size_t> const &nc_start,    // Where to start each dimension in NetCDF
    std::vector<int> const &b2n,    // Where to slot each Blitz++ dimension
    std::vector<std::string> const &arr_sdims)
{
    netCDF::NcVar ncvar = ncio.nc->getVar(vname);
    // By necessity, we must list dimensions in NetCDF order, so we know
    // where Blitz++ dimensions slot in
    bool const ncdims_in_nc_order=true;

    auto ordering(tiny_to_vector(arr.ordering()));
    auto ncdims_nd(named_dims(ncdims, ncdims_in_nc_order, ordering, ncvar));

    _ncio_blitz::Info info(
        named_dims(arr, arr_sdims), ordering,
        DimOrderMatch::MEMORY,    // Not used
        std::move(ncdims_nd),
        std::vector<int>(b2n), std::vector<size_t>(nc_start));

    return _ncio_blitz::ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS, info);
}

#undef NCIO_BLITZ_PARAMS
#undef NCIO_BLITZ_ARGS

// =================================================

/** Convenience method: reads data from NetCDF, returns as a newly allocate blitz::Array */
template<class TypeT, int RANK>
inline blitz::Array<TypeT,RANK> nc_read_blitz(
    netCDF::NcGroup *nc,
    std::string const &vname,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>())
{
    NcIO ncio(nc, 'r');    // Dummy
    blitz::Array<TypeT,RANK> val;
    ncio_blitz_alloc(ncio, val, vname, get_nc_type<TypeT>(), {}, DimOrderMatch::MEMORY, true, storage);
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

    return std::vector<netCDF::NcDim>
        { get_or_add_dim(ncio, dim_names[0], val.size()) };
}



template<class TypeT, class CastT>
void nc_rw_vector(
    netCDF::NcGroup *nc,
    char rw,
    std::vector<TypeT> *val,
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
            ncvar.getVar(startp, countp, (CastT *)&(*val)[0]);
        break;
        case 'w' :
            ncvar.putVar(startp, countp, (CastT *)&(*val)[0]);
        break;
    }
}

template<class TypeT, class CastT=TypeT>
std::vector<TypeT> nc_read_vector(
    netCDF::NcGroup *nc,
    std::string const &vname)
{
    std::vector<TypeT> val;
    nc_rw_vector(nc, 'r', &val, true, vname);
    return val;
}


/** Define and write a std::vector. */
template<class TypeT, class CastT=TypeT>
netCDF::NcVar ncio_vector(
    NcIO &ncio,
    std::vector<TypeT> &val,
    bool alloc,         // Should we allocate val?
    std::string const &vname,
    std::string const &snc_type,
    std::vector<netCDF::NcDim> const &dims)
{
    netCDF::NcVar ncvar = get_or_add_var(ncio, vname, snc_type, dims);

    ncio += std::bind(&nc_rw_vector<TypeT,CastT>, ncio.nc, ncio.rw, &val, alloc, vname);
    return ncvar;
}
// ----------------------------------------------------

// ----------------------------------------------------

/** Do linewrap for strings that are intended to be used as comment attributes in NetCDF files.
       see: http://www.cplusplus.com/forum/beginner/19034/
*/
extern std::string ncwrap( std::string const &str, size_t width = 55 );


// ===========================================================
/** Stores the value of any NetCDF attribute, of any type.
This is useful for copying attributes from one variable to another.
We don't yet have a way to extract regular values from these. */
class NcAttValue {
    std::unique_ptr<netCDF::NcType> nctype;
    std::vector<char> data;
public:
    size_t getAttLength() const
        { return data.size() / nctype->getSize(); }
    size_t getTypeSize() const
        { return nctype->getSize(); }

    static NcAttValue get(netCDF::NcAtt const &ncatt);
    void put(netCDF::NcVar const &ncvar, std::string const &name) const;
};


extern void get_or_put_att(
    netCDF::NcVar &ncvar, char rw,
    const std::string &name,
    NcAttValue &aval);


extern netCDF::NcVar get_or_put_all_atts(
    netCDF::NcVar &ncvar, char rw,
    std::vector<std::pair<std::string, NcAttValue>> &avals);

extern std::vector<std::pair<std::string, NcAttValue>> get_all_atts(netCDF::NcVar const &ncvar);

}   // Namespace
#endif  // Guard
