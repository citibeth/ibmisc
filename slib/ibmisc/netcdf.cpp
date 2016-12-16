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
#include <ibmisc/netcdf.hpp>


using namespace netCDF;

namespace ibmisc {

bool netcdf_debug = false;


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

NcIO::NcIO(std::string const &filePath, char mode) :
    _mync(filePath, _filemode_to_netcdf(mode),
        netCDF::NcFile::FileFormat::nc4),
    own_nc(true),
    nc(&_mync),
    rw(_filemode_to_rw(mode)),
    define(rw == 'w') {}

void NcIO::operator+=(std::function<void ()> const &fn)
{
    if (rw == 'r') fn();
    else _io.push_back(fn);
}

void NcIO::operator()() {
    for (auto ii=_io.begin(); ii != _io.end(); ++ii) (*ii)();
    _io.clear();
    tmp.free();
}

void NcIO::close() {
    if (own_nc) {
        (*this)();
        _mync.close();    // NcFile::close() is idempotent
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

    if (ncio.rw == 'w' && dim.getSize() != dim_size) {
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
        std::vector<std::string> sdims;
        for (auto dim=dims.begin(); dim != dims.end(); ++dim)
            sdims.push_back(dim->getName());

        if (ncvar.isNull()) {
            ncvar = ncio.nc->addVar(vname, snc_type, sdims);
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

}   // Namespace

