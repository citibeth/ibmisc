#include <ibmisc/ncfile.hpp>
#include <ibmisc/blitz.hpp>

using namespace netCDF;

namespace ibmisc {

#if 0
static void nc_write_time0(netCDF::NcGroup *nc, double time0, std::string const &vname)
{
    netCDF::NcVar time0_var(nc->getVar(vname));
    time0_var.putVar({0}, {1}, &time0);
}

void ncdefine_time0(NcIO &ncio, double time0, std::string const &time_units, std::string const &vname_base)
{
    // if (ncio.rw == 'r') (*ibmisc_error)(-1, "Read not supported");

    NcDim time_dim(get_or_add_dim(ncio, "time", -1));

    NcVar time0_var = ncio.nc->addVar(vname_base + "time0", ibmisc::get_nc_type<double>(), dims_b);
    time0_var.putAtt("units", time_units);
    time0_var.putAtt("calendar", "365_day");
    time0_var.putAtt("axis", "T");
    time0_var.putAtt("long_name", "Simulation start time");

    NcVar time_var = ncio.nc->addVar(vname_base + "time", ibmisc::get_nc_type<double>(), dims_atmosphere);
    time_var.putAtt("units", time_units);
    time_var.putAtt("calendar", "365_day");
    time_var.putAtt("axis", "T");
    time_var.putAtt("long_name", "Coupling times");

    ncio += std::bind(&nc_write_time0, ncio.nc, time0_s, vname_base + "time0");
}
#endif
// ----------------------------------------------------------
static void nc_write_timespan(netCDF::NcGroup *nc, std::array<double,2> const &timespan, std::string const &vname)
{
    netCDF::NcVar timespan_var(nc->getVar(vname));
    timespan_var.putVar({0}, {2}, &timespan[0]);
}

/** For files with a single timepoint: Defines start and end timepoints. */
void ncio_timespan(NcIO &ncio, std::array<double,2> &timespan, std::string const &time_units, std::string const &vname_base)
{
    auto timespan_b(to_blitz<double,2>(timespan));
    std::string vname = vname_base + "timespan";
    ncio_blitz<double,1>(ncio, timespan_b, false, vname, "double",
        get_or_add_dims(ncio, {"two"}, {2}));

    if (ncio.rw == 'w') {
        NcVar timespan_var(ncio.nc->getVar(vname));
        timespan_var.putAtt("units", time_units);
        timespan_var.putAtt("calendar", "365_day");
        timespan_var.putAtt("axis", "T");
        timespan_var.putAtt("long_name", "[start, now] of this time segment");
    }
}
// ----------------------------------------------------------

}
