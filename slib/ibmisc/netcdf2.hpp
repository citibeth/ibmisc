#ifndef IBMISC_NETCDF2_HPP
#define IBMISC_NETCDF2_HPP

#include <ibmisc/netcdf.hpp>
#include <functional>

/** Alternative ways to read/write Blitz++ in NetCDF */

namespace ibmisc {

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
// -------------------------------------------------------------------
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
template<int RANK>
struct NcIOBlitzInfo {
    std::vector<netCDF::NcDim> dims;    // Dimensions of ON-DISK NetCDF Variable
    std::vector<size_t> nc_start;        // Where to the NetCDF variable; same rank as dims.   If {}, then all 0's.  <0 means "allocate a Blitz var for this dim."
    std::array<int,RANK> b2n;    // (Blitz dim i) corresponds to (NetCDF dim b2n[i]).  If {}, then {0,1,2,...}; otherise, rank=RANK
};

// ------------------------------------------------------------------------------
template<class TypeT, int RANK>
using ncio_blitz_helper = std::function<NcIOBlitzInfo<RANK> (NcIO &, netCDF::NcVar &, blitz::Array<TypeT,RANK> &)>;

#define NCIO_BLITZ_PARAMS \
    NcIO &ncio, \
    blitz::Array<TypeT, RANK> &arr, \
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
    ncio_blitz_helper<TypeT,RANK> const &info_fn);

/** Maps a blitz::Array directly to/from a NetCDF array */
template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz(
    NCIO_BLITZ_PARAMS,
    ncio_blitz_helper<TypeT,RANK> const &info_fn)
{
    // nc->getVar() returns a null NcVar if no object of that name is found.
    // That is what we want here.
    netCDF::NcVar ncvar(ncio.nc->getVar(vname));
    NcIOBlitzInfo<RANK> const info(info_fn(ncio, ncvar, arr));

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
            "User-supplied dimension[%d] does not match netCDF dimension=%s",
            k, ncvar.getDim(k).getName().c_str());
    }

    // Check: in-memory variable is compatible with NetCDF variable
    for (int kb=0; kb<RANK; ++kb) {
        int const kn = info.b2n[kb];
        if (arr.extent(kb) > ncvar.getDim(kn).getSize() - info.nc_start[kn]) (*ibmisc_error)(-1,
            "C++ variable extent[%d]=%d does not match NetCDF variable extent[%d]=%d",
            kb, arr.extent(kb), kn, ncvar.getDim(info.b2n[kn]).getSize(), info.nc_start[kn]);
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
NcIOBlitzInfo<RANK> _whole1(
    // Leave these unbound
    NcIO &ncio,
    netCDF::NcVar &ncvar,
    blitz::Array<TypeT, RANK> &arr,
    // Bind these away
    std::vector<netCDF::NcDim> _dims,
    bool equal_dim_order,
    bool dims_in_nc_order);

template<class TypeT, int RANK>
NcIOBlitzInfo<RANK> _whole1(
    // Leave these unbound
    NcIO &ncio,
    netCDF::NcVar &ncvar,
    blitz::Array<TypeT, RANK> &arr,
    // Bind these away
    std::vector<netCDF::NcDim> _dims,
    bool equal_dim_order,
    bool dims_in_nc_order)
{
    int const nc_rank = RANK;

    NcIOBlitzInfo<RANK> ret;

    // If no dimensions provided, we need to get them from NetCDF array
    if (_dims.size() == 0) _dims = ncvar.getDims();

    // Check ranks
    if (_dims.size() != RANK) (*ibmisc_error)(-1,
        "Rank of NetCDF dims (=%d) must correspond to Blitz rank (=%d)",
        (int)_dims.size(), RANK);

    // Put ret.dims in NetCDF order
    if (!dims_in_nc_order && !equal_dim_order) {
        // Re-order user-supplied dimensions to NetCDF order
        // (that is, decreasing rank)
        for (int in=0; in<nc_rank; ++in) {
            int const ib = arr.ordering(RANK-in-1);
            ret.dims.push_back(_dims[ib]);
        }
    } else {
        // Dimensions already in NetCDF order, just copy them.
        ret.dims = std::move(_dims);
    }

    // Set ret.b2n
    if (equal_dim_order) {
        for (int ib=0; ib<RANK; ++ib) ret.b2n[ib] = ib;
    } else {
        for (int in=0; in<nc_rank; ++in) {
            int const ib = arr.ordering(RANK-in-1);
            ret.b2n[ib] = in;
        }
    }

    // Zero out start, we're doing whole array here...
    for (int i=0; i<nc_rank; ++i) ret.nc_start.push_back(0);

    return ret;
}

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_whole(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &dims,
    bool equal_dim_order=false,
    bool dims_in_nc_order=true)
{
    using namespace std::placeholders;
    return ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS,
        std::bind(&_whole1<TypeT,RANK>, _1, _2, _3,
            dims, equal_dim_order, dims_in_nc_order));
}
// -------------------------------------------------------------------------------
template<class TypeT, int RANK>
NcIOBlitzInfo<RANK> _alloc1(
    // Leave these unbound
    NcIO &ncio,
    netCDF::NcVar &ncvar,
    blitz::Array<TypeT, RANK> &arr,
    // Bind these away
    std::vector<netCDF::NcDim> _dims,
    bool dims_in_nc_order,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    // Array in memory must (will be allocated to) have same
    // memory layout as on disk.
    bool const equal_dim_order = false;

    // Determine dimensions
    if (_dims.size() == 0) {
        _dims = ncvar.getDims();
        dims_in_nc_order = true;
    }

    // Allocate the Blitz array
    // This only makes sense for reading
    if (ncio.rw == 'r') {

        blitz::TinyVector<int,RANK> extent;
        for (int in=0; in<RANK; ++in) {
            int ib = storage.ordering(RANK-in-1);
            extent[ib] = _dims[in].getSize();
        }
        arr.reference(blitz::Array<TypeT,RANK>(extent, storage));
    }

    return _whole1(ncio, ncvar, arr,
        _dims, equal_dim_order, dims_in_nc_order);
}

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_alloc(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> _dims = {},
    bool dims_in_nc_order = true,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>())
{
    using namespace std::placeholders;
    return ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS,
        std::bind(&_alloc1<TypeT,RANK>, _1, _2, _3,
            _dims, dims_in_nc_order, storage));
}
// -------------------------------------------------------------------------
/** Convenience method: reads data from NetCDF, returns as a newly allocate blitz::Array */
template<class TypeT, int RANK>
blitz::Array<TypeT,RANK> nc_read_blitz(
    netCDF::NcGroup *nc,
    std::string const &vname,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

template<class TypeT, int RANK>
blitz::Array<TypeT,RANK> nc_read_blitz(
    netCDF::NcGroup *nc,
    std::string const &vname,
    blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>())
{
    NcIO ncio(nc, 'r');    // Dummy
    blitz::Array<TypeT,RANK> val;
    ncio_blitz_alloc(ncio, val, vname, get_nc_type<TypeT>(), {}, true, storage);
    return val;
}
// -------------------------------------------------------------------------
template<class TypeT, int RANK>
NcIOBlitzInfo<RANK> _partial(
    // Leave these unbound
    NcIO &ncio,
    netCDF::NcVar &ncvar,
    blitz::Array<TypeT, RANK> &arr,
    // Bind these away
    std::vector<netCDF::NcDim> dims,
    std::vector<int> nc_start,    // Where to start each dimension in NetCDF
    std::array<int,RANK> const &b2n)    // Where to slot each Blitz++ dimension
;

template<class TypeT, int RANK>
NcIOBlitzInfo<RANK> _partial(
    // Leave these unbound
    NcIO &ncio,
    netCDF::NcVar &ncvar,
    blitz::Array<TypeT, RANK> &arr,
    // Bind these away
    std::vector<netCDF::NcDim> dims,
    std::vector<int> nc_start,    // Where to start each dimension in NetCDF
    std::array<int,RANK> const &b2n)    // Where to slot each Blitz++ dimension
{
    // By necessity, we must list dimensions in NetCDF order, so we know
    // where Blitz++ dimensions slot in
    bool const dims_in_nc_order=true;

    NcIOBlitzInfo<RANK> ret;

    // Check ranks
    if (dims.size() < RANK) (*ibmisc_error)(-1,
        "Rank of NetCDF dims (=%d) must be at least Blitz rank (=%d)",
        (int)dims.size(), RANK);

    // Put ret.dims in NetCDF order
    ret.dims = std::move(dims);

    // Set ret.nc_start
    ret.nc_start = std::move(nc_start);

    // Set ret.b2n
    ret.b2n = b2n;

    return ret;
}

template<class TypeT, int RANK>
netCDF::NcVar ncio_blitz_partial(
    NCIO_BLITZ_PARAMS,
    std::vector<netCDF::NcDim> const &dims,
    std::vector<int> const &nc_start,    // Where to start each dimension in NetCDF
    std::array<int,RANK> const &b2n)    // Where to slot each Blitz++ dimension
{
    using namespace std::placeholders;
    return ncio_blitz<TypeT,RANK>(NCIO_BLITZ_ARGS,
        std::bind(&_partial<TypeT,RANK>, _1, _2, _3,
            dims, nc_start, b2n));
}




}
#endif
