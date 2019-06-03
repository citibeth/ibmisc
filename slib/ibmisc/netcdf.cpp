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

#include <string>
#include <netcdf>
#include <sstream>
#include <boost/filesystem.hpp>
#include <ibmisc/netcdf.hpp>


using namespace netCDF;

namespace ibmisc {

bool netcdf_debug = false;

netCDF::NcType nc_type(netCDF::NcVar const &ncvar, std::string sntype)
{
    netCDF::NcType ret(ncvar.getParentGroup().getType(sntype, netCDF::NcGroup::ParentsAndCurrent));

    // See legal data types:
    //     https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html
    if (ret.isNull()) (*ibmisc_error)(-1,
        "Illegal NetCDF type \"%s\".  Valid types are probably byte, ubyte, char, short, ushort, int, uint, int64, uint64, float, double, string", sntype.c_str());
    return ret;
}

netCDF::NcType nc_type(netCDF::NcGroup const &ncgroup, std::string sntype)
{
    netCDF::NcType ret(ncgroup.getType(sntype, netCDF::NcGroup::ParentsAndCurrent));

    // See legal data types:
    //     https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html
    if (ret.isNull()) (*ibmisc_error)(-1,
        "Illegal NetCDF type \"%s\".  Valid types are probably byte, ubyte, char, short, ushort, int, uint, int64, uint64, float, double, string", sntype.c_str());
    return ret;
}



void _check_nc_rank(
    netCDF::NcVar const &ncvar,
    int rank)
{
    // Check rank of NetCDF variable
    if (ncvar.getDimCount() != rank)
        (*ibmisc_error)(-1,
            "NetCDF variable (%s) of rank %ld does not match blitz::Array "
            "of rank %d", ncvar.getName().c_str(), ncvar.getDimCount(), rank);
}

std::vector<NamedDim> named_dims(netCDF::NcVar &ncvar)
{
    std::vector<NamedDim> ret;
    if (ncvar.isNull()) return ret;

    for (int in=0; in<ncvar.getDimCount(); ++in) {
        netCDF::NcDim ncdim(ncvar.getDim(in));
        ret.push_back(NamedDim(ncdim.getName(), ncdim.getSize()));
    }
    return ret;
}

std::vector<NamedDim> named_dims(
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<int> const &ordering)
{
    bool const nc_dims_in_nc_order = (ordering.size() == 0);
    std::vector<NamedDim> ret;

    for (int in=0; in<ncdims.size(); ++in) {
        int const ib = (nc_dims_in_nc_order ? in : ordering[ordering.size()-in-1]);
//printf("in=%d, ib=%d\n", in, ib);
        netCDF::NcDim const &ncdim(ncdims[ib]);

        if (ncdim.isNull()) {
            ret.push_back(NamedDim());    // slot / placeholder
        } else {
            ret.push_back(NamedDim(ncdim.getName(), ncdim.getSize()));
        }
    }
    return ret;
}

std::vector<NamedDim> named_dims(
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<int> const &ordering,
    netCDF::NcVar &ncvar)
{
    if (ncdims.size() == 0) {
        auto ret(named_dims(ncvar));
#if 0
        // This is just not true.  See Info::Info(), named_dims result can have zero size.
        // This check maybe improves error messages in some places. But it also breaks other stuff.
        // Figure out a better way to do it.
        if (ret.size() == 0) (*ibmisc_error)(-1,
            "named_dims() found nothing.  NOTE: ncio_blitz() writes require dimensions to be supplied.");
#endif
        return ret;
    } else {
        return named_dims(ncdims, ordering);
    }
}

std::vector<NamedDim> named_dims(
    std::vector<netCDF::NcDim> const &dims,
    bool dims_in_nc_order,
    std::vector<int> const &blitz_ordering,
    netCDF::NcVar &ncvar)
{
    if (dims_in_nc_order) {
        return named_dims(dims, std::vector<int>{}, ncvar);
    } else {
        return named_dims(dims, blitz_ordering, ncvar);
    }
}


/** Reconciles dimensions between Blitz and NetCDF.  This subroutine
    fills in the blanks, it does not check afterwards.  That is for elsewhere. */
_ncio_blitz::Info::Info(
    std::vector<NamedDim> &&_blitz,
    std::vector<int> const &blitz_ordering,    // Ordering of dimensions in increasing stride
    DimOrderMatch match,
    std::vector<NamedDim> &&_netcdf,
    std::vector<int> &&_b2n,
    std::vector<size_t> &&_nc_start)
: blitz(std::move(_blitz)), netcdf(_netcdf), b2n(std::move(_b2n)), nc_start(std::move(_nc_start))
{
    int const the_case =
        (blitz.size() == 0 ? 0 : 4) +
        (netcdf.size() == 0 ? 0 : 2) +
        (b2n.size() == 0 ? 0 : 1);
//printf("the_case = %d\n", the_case);
    switch(the_case) {
        case 0:    // NOTHING
        case 1: {
                (*ibmisc_error)(-1,
                    "Unable to reconcile blitz + netcdf when nothing is given");
        } break;
        case 2: {    // netcdf (blitz_ordering, match)
            int const netcdf_rank = netcdf.size();
            int const blitz_rank = netcdf_rank;
            for (int ib=0; ib<blitz_rank; ++ib) {
                b2n.push_back(-1);
                blitz.push_back(NamedDim());
            }
            for (int in=0; in<netcdf_rank; ++in) {
                int const ib = (match == DimOrderMatch::LEXICAL ? in : blitz_ordering[blitz_rank-in-1]);
                blitz[ib] = netcdf[in];
                b2n[ib] = in;
            }
        } break;
        case 3: {    // netcdf, b2n
            int const blitz_rank = b2n.size();
            for (int ib=0; ib<blitz_rank; ++ib) {
                int const in = b2n[ib];
                blitz.push_back(netcdf[in]);
            }
        } break;
        case 4: {    // blitz (blitz_ordering, match)
            int const blitz_rank = blitz.size();
            int const netcdf_rank = blitz_rank;

            b2n.clear();
            for (int ib=0; ib<blitz_rank; ++ib) b2n.push_back(-1);
            for (int in=0; in<netcdf_rank; ++in) {
                int const ib = (match == DimOrderMatch::LEXICAL ? in : blitz_ordering[blitz_rank-in-1]);
                b2n[ib] = in;
                netcdf.push_back(blitz[ib]);
            }
        } break;
        case 5: {    // blitz, b2n
            for (int in=0; in<blitz.size(); ++in) netcdf.push_back(NamedDim());

            for (int ib=0; ib<blitz.size(); ++ib) {
                int const in = b2n[ib];
                netcdf[in] = blitz[ib];
            }
        } break;
        case 6: {    // blitz, netcdf, (blitz_ordering, match)
            int const blitz_rank = blitz.size();
            int const netcdf_rank = netcdf.size();

            int nname_blitz = 0;
            for (int ib=0; ib<blitz_rank; ++ib)
                if (blitz[ib].name != "") ++nname_blitz;
            if (nname_blitz != 0 && nname_blitz != blitz_rank)
                (*ibmisc_error)(-1, "Blitz dimensions must all have names or none have names");

            int nname_netcdf = 0;
            for (int in=0; in<netcdf_rank; ++in) {
                if (netcdf[in].name != "") ++nname_netcdf;
            }
            if (nname_netcdf != netcdf_rank - blitz_rank && nname_netcdf != netcdf_rank)
                (*ibmisc_error)(-1, "NetCDF dimensions must all have names or none have names");

            int name_case =
                (nname_blitz == blitz_rank ? 2 : 0) +
                (nname_netcdf == netcdf_rank ? 1 : 0);
//printf("name_case = %d\n", name_case);
            switch(name_case) {
                case 0: (*ibmisc_error)(-1,
                    "Names must be specified for at least Blitz or NetCDF variable");
                break;
                case 3 : {
                    // Names are fully specified: match by names, ignore slots

                    // Set up name lookup
                    std::map<std::string, int> name2n;
                    for (int in=0; in<netcdf_rank; ++in) {
                        std::string const &name(netcdf[in].name);
                        if (name.size() > 0)
                            name2n.insert(std::make_pair(netcdf[in].name, in));
                    }
                    for (int ib=0; ib<blitz_rank; ++ib) {
                        try {
                            int const in = name2n.at(blitz[ib].name);
                            b2n.push_back(in);
                        } catch(std::out_of_range const &e) {
                            (*ibmisc_error)(-1, "_ncio_blitz::Info(): Error looking up dimension name %s\n", blitz[ib].name.c_str());
                        }
                    }
                } break;
                case 1:        // Names only for NetCDF
                case 2 : {   // Names only for blitz.
                    // Names only specified for one of Blitz or NetCDF

                    // ------- Set b2n
                    for (int ib=0; ib<blitz_rank; ++ib) b2n.push_back(-1);

                    // Count slots in NetCDF
                    int nslot = 0;
                    for (int in_full=0; in_full<netcdf_rank; ++in_full)
                        if (!netcdf[in_full].extent < 0) ++nslot;

                    if (nslot == 0) {
                        // No slots in NetCDF; match dimensions by
                        // ordering and hope for the best.
                        for (int in=0; in<netcdf_rank; ++in) {
                            int const ib = (match == DimOrderMatch::LEXICAL ? in
                                : blitz_ordering[blitz_rank-in-1]);
                            b2n[ib] = in;
                        }
                    } else if (nslot == blitz_rank) {
                        // Fill Blitz dimensions into slots in NetCDF
                        int in_sub = 0;
                        for (int in_full=0; in_full<netcdf_rank; ++in_full) {
                            if (!netcdf[in_full].extent < 0) {
                                int const ib = (match == DimOrderMatch::LEXICAL ? in_sub
                                    : blitz_ordering[blitz_rank-in_sub-1]);
                                b2n[ib] = in_full;
                                ++in_sub;
                            }
                        }
                    } else {
                        (*ibmisc_error)(-1,
                            "Number of slots=%d must equal Blitz rank=%d, or zero", nslot, blitz_rank);
                    }
                } break;
            }    // switch(name_case)
        }    // case 6: blitz, netcdf specified
        case 7: {    // blitz, netcdf, b2n
            // Fully specified; nothing to do here
        } break;
    }    // switch(the_case)

    // Fix names and extents, as necessary
    for (int ib=0; ib<b2n.size(); ++ib) {
        int const in = b2n[ib];

        if (in < 0 || in >= netcdf.size()) (*ibmisc_error)(-1,
            "At Blitz index=%d, NetCDF inex=%d is out of range [0,%ld)", ib, in, netcdf.size());

        if (netcdf[in].name == "") {
            netcdf[in].name = blitz[ib].name;
        } else if (blitz[ib].name == "") {
            blitz[ib].name = netcdf[in].name;
        } else if (blitz[ib].name != netcdf[in].name) {
            (*ibmisc_error)(-1,
                "Dimension name mismatch: blitz[%d]=%s, netcdf[%d]=%s",
                ib, blitz[ib].name.c_str(),
                in, netcdf[in].name.c_str());
        }

        if (netcdf[in].extent < 0) {
            netcdf[in].extent = blitz[ib].extent;
        } else if (blitz[ib].extent < 0) {
            blitz[ib].extent = netcdf[in].extent;
        } else if (blitz[ib].extent != netcdf[in].extent) {
            (*ibmisc_error)(-1,
                "Dimension extent mismatch: blitz[dim=%d:%s]=%d, netcdf[dim=%d:%s]=%d",
                ib, blitz[ib].name.c_str(), blitz[ib].extent,
                in, netcdf[in].name.c_str(), netcdf[in].extent);
        }
    }

    if (nc_start.size() == 0)
        for (int i=0; i<netcdf.size(); ++i)
            nc_start.push_back(0);
}

// ===============================================

static std::map<std::string, netCDF::NcFile::FileFormat> const _sformat_to_netcdf {
    std::make_pair("classic", NcFile::FileFormat::classic),
    std::make_pair("classic64", NcFile::FileFormat::classic64),
    std::make_pair("nc4", NcFile::FileFormat::nc4),
    std::make_pair("nc4classic", NcFile::FileFormat::nc4classic)
};

/*
@param mode Can be (see https://docs.python.org/3/library/functions.html#open)
 'r' open for reading (default)
 'w' open for writing, truncating the file first
 'x' open for exclusive creation, failing if the file already exists
 'a' open for writing, appending to the end of the file if it exists
*/
inline netCDF::NcFile::FileMode _filemode_to_netcdf(char mode)
{
    switch(mode) {
        case 'r' : return netCDF::NcFile::FileMode::read;
        case 'w' : return netCDF::NcFile::FileMode::replace;
        case 'x' : return netCDF::NcFile::FileMode::newFile;
        case 'a' : return netCDF::NcFile::FileMode::write;

        // Hack: the netCDF labels can also be coerced to char
        case netCDF::NcFile::FileMode::read:
        case netCDF::NcFile::FileMode::replace:
        case netCDF::NcFile::FileMode::newFile:
        case netCDF::NcFile::FileMode::write:
            return (netCDF::NcFile::FileMode)mode;
    }
    (*ibmisc_error)(-1,
        "Illegal filemode: '%c'", mode);
}
inline char _filemode_to_rw(char mode)
{
    switch(mode) {
        case 'r' :
        case netCDF::NcFile::FileMode::read:
            return 'r';
        default:
            return 'w';
    }
}

void NcIO::default_configure_var(netCDF::NcVar ncvar)
{
    ncvar.setCompression(true, true, 4);

    // For some reason, this causes an HDF5 error
    // ncvar.setChecksum(netCDF::NcVar::nc_FLETCHER32);
}

void NcIO::no_compress(netCDF::NcVar ncvar)
{
    ncvar.setCompression(true, true, 4);

    // For some reason, this causes an HDF5 error
    // ncvar.setChecksum(netCDF::NcVar::nc_FLETCHER32);
}

inline std::unique_ptr<NcFile> open_netcdf(std::string const &filePath, char mode, std::string const &sformat)
{
    NcFile::FileMode fmode(_filemode_to_netcdf(mode));
    switch(fmode) {
        case NcFile::FileMode::read:
        case NcFile::FileMode::write:
            if (!boost::filesystem::exists(filePath)) (*ibmisc_error)(-1,
                "Trying to open file %s, which doesn't exist", filePath.c_str());
        break;
    }


    netCDF::NcFile::FileFormat format = _sformat_to_netcdf.at(sformat);
    return std::unique_ptr<NcFile>(new NcFile(filePath, fmode, format));
}

NcIO::NcIO(std::string const &filePath, char mode,
    std::string const &sformat,
    std::function<void(NcVar)> const &_configure_var) :
    _mync(open_netcdf(filePath, mode, sformat)),
    nc(&*_mync),
    rw(_filemode_to_rw(mode)),
    define(rw == 'w'),
    configure_var(_configure_var)
{
}

void NcIO::add(std::string const &tag, std::function<void ()> const &fn)
{
    if (rw == 'r') fn();
    else _io.push_back(TaggedThunk(fn, tag));
}

void NcIO::flush(bool debug) {
    for (auto ii=_io.begin(); ii != _io.end(); ++ii) {
        if (debug) printf("NcIO::flush(%s)\n", ii->tag.c_str());
        ii->fn();
    }
    _io.clear();
    tmp.free();
}

void NcIO::close() {
    if (_mync.get()) {
        (*this)();
        _mync.reset();
    }
}


// -----------------------------------------------------
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
            "Attempt in get_or_add_dim(%s) to change size from %d to unlimited",
            dim_name.c_str(), dim.getSize());
    }

    return dim;
}

netCDF::NcDim get_or_add_dim(NcIO &ncio, std::string const &dim_name, size_t dim_size)
{
    // dim_size < 0 means "unlimited"
    if (dim_size < 0) return get_or_add_dim(ncio, dim_name);

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

    if (ncio.rw == 'w' && dim.getSize() != dim_size) {
        (*ibmisc_error)(-1, 
            "Attempt in get_or_add_dim(%s) to change size from %ld to %ld",
            dim_name.c_str(), dim.getSize(), dim_size);
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

std::vector<netCDF::NcDim> get_or_add_finite_dims(
    NcIO &ncio,
    std::vector<std::string> const &dim_names,
    std::vector<size_t> const &dim_sizes)
{
    std::vector<long> sizes;
    sizes.reserve(dim_sizes.size());
    for (size_t s : dim_sizes) sizes.push_back((long)s);
    return get_or_add_dims(ncio, dim_names, sizes);
}

std::vector<netCDF::NcDim> get_or_add_dims(
    NcIO &ncio,
    std::vector<std::string> const &dim_names,
    std::vector<long> const &dim_sizes)
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
            // Unlimited
            ret[k] = get_or_add_dim(ncio, dim_names[k]);
        } else {
            // Regular
            ret[k] = get_or_add_dim(ncio, dim_names[k], dim_sizes[k]);
        }
    }
    return ret;
}
// ====================================================
/** This works for all valid NetCDF types. */
netCDF::NcVar get_or_add_var(
    NcIO &ncio,
    std::string const &vname,
    std::string const &snc_type,
    std::vector<netCDF::NcDim> const &dims)
{
    netCDF::NcVar ncvar;
    if (ncio.define) {
        ncvar = ncio.nc->getVar(vname);

        if (ncvar.isNull()) {
            NcType nctype(nc_type(*ncio.nc, snc_type));
            ncvar = ncio.nc->addVar(vname, nctype, dims);
            ncio.configure_var(ncvar);
        } else {
            // Check dimensions match
            if (ncvar.getDimCount() != dims.size()) {
                (*ibmisc_error)(-1,
                    "NetCDF variable %s(%d dims) has wrong number of "
                    "dimensions, %d expected",
                    vname.c_str(), ncvar.getDimCount(), dims.size());
            }
            for (int i=0; i<ncvar.getDimCount(); ++i) {
                NcDim ncdim = ncvar.getDim(i);
                if (ncdim != dims[i]) {
                    (*ibmisc_error)(-1,
                        "Trying to change dimension %d of "
                        "NetCDF variable %s from %s=%ld to %s=%ld",
                        i, ncvar.getName().c_str(),
                        ncdim.getName().c_str(), ncdim.getSize(),
                        dims[i].getName().c_str(), dims[i].getSize());
                }
            }
        }
    } else {
        ncvar = ncio.nc->getVar(vname);
        if (ncvar.isNull()) {
            (*ibmisc_error)(-1,
                "Variable %s required but not found", vname.c_str());
        }

        // Check that types match
        NcType spec_type(nc_type(ncvar, snc_type));
        NcType real_type(ncvar.getType());
        if (spec_type != real_type) (*ibmisc_error)(-1,
            "On-disk type of variable %s = %s does not match specified type %s (%s)",
            vname.c_str(), real_type.getName().c_str(),
            snc_type.c_str(), spec_type.getName().c_str());
    }
    return ncvar;
}

// ---------------------------------------------
netCDF::NcVarAtt get_att(netCDF::NcVar &ncvar, std::string const &name)
{
    try {
        return ncvar.getAtt(name);
    } catch(netCDF::exceptions::NcException &e) {
        (*ibmisc_error)(-1,
            "Attribute %s:%s does not exist",
            ncvar.getName().c_str(), name.c_str());
    }
}


netCDF::NcGroupAtt get_att(netCDF::NcGroup &ncvar, std::string const &name)
{
    try {
        return ncvar.getAtt(name);
    } catch(netCDF::exceptions::NcException &e) {
        (*ibmisc_error)(-1,
            "Attribute %s:%s does not exist",
            ncvar.getName().c_str(), name.c_str());
    }
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

// ==============================================

NcAttValue NcAttValue::get(NcAtt const &ncatt)
{
    NcAttValue aval;
    aval.nctype.reset(new NcType(ncatt.getType()));
    size_t const nbytes = ncatt.getAttLength() * aval.nctype->getSize();
    aval.data.resize(nbytes);
    ncatt.getValues((void *)&aval.data[0]);
    return aval;
}
void NcAttValue::put(NcVar const &ncvar, std::string const &name) const
{
    ncvar.putAtt(name, *nctype, getAttLength(), &data[0]);
}


void get_or_put_att(
    NcVar &ncvar, char rw,
    const std::string &name,
    NcAttValue &aval)
{
    if (rw == 'r') {
        NcVarAtt ncatt(ncvar.getAtt(name));
        aval = NcAttValue::get(ncatt);
    } else {
        aval.put(ncvar, name);
    }
}


NcVar get_or_put_all_atts(
    NcVar &ncvar, char rw,
    std::vector<std::pair<std::string, NcAttValue>> &avals)
{
    if (rw == 'r') {
        avals.clear();
        auto atts(ncvar.getAtts());
        for (auto ii=atts.begin(); ii != atts.end(); ++ii) {
            avals.push_back(std::make_pair(
                ii->first, NcAttValue::get(ii->second)));

        }
    } else {
        for (auto ii=avals.begin(); ii != avals.end(); ++ii) {
            std::string const &name(ii->first);
            NcAttValue const &aval(ii->second);

            aval.put(ncvar, name);
        }
    }
    return ncvar;
}

std::vector<std::pair<std::string, NcAttValue>> get_all_atts(NcVar const &ncvar)
{
    std::vector<std::pair<std::string, NcAttValue>> ret;
    get_or_put_all_atts(*const_cast<NcVar *>(&ncvar), 'r', ret);
    return ret;
}



}   // Namespace

