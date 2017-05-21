#include <ibmisc/datetime.hpp>
#include <ibmisc/ncfile.hpp>
#include <ibmisc/blitz.hpp>

using namespace netCDF;

namespace ibmisc {

// ----------------------------------------------------------
static void nc_write_timespan(netCDF::NcGroup *nc, std::array<double,2> const &timespan, std::string const &vname)
{
    netCDF::NcVar timespan_var(nc->getVar(vname));
    timespan_var.putVar({0}, {2}, &timespan[0]);
}

/** For files with a single timepoint: Defines start and end timepoints. */
void ncio_timespan(NcIO &ncio, std::array<double,2> &timespan, ibmisc::TimeUnit const &time_unit, std::string const &vname)
{
    auto &timespan_b(ncio.tmp.move(to_blitz<double,2>(timespan)));    // Copy the input
    ncio_blitz<double,1>(ncio, timespan_b, false, vname, "double",
        get_or_add_dims(ncio, {vname + ".length"}, {timespan.size()}));

    if (ncio.rw == 'w') {
        // Write human-readable form of the timestamps
        auto &timespan_c(ncio.tmp.move(    // Allocate for the NcIO operation
            blitz::Array<char,2>(timespan.size(), iso8601_length)));

        for (int i=0; i < timespan.size(); ++i) {
            ibmisc::Datetime dt(time_unit.to_datetime(timespan[i]));
            std::string sdate = to_iso8601(dt);
            for (int j=0; j<iso8601_length; ++j) timespan_c(i,j) = sdate[j];
        }
        ncio_blitz<char,2>(ncio, timespan_c, false, vname+".txt", "char",
            get_or_add_dims(ncio,
                {vname + ".length", "iso8601.length"},
                {timespan.size(), iso8601_length}));

        // Add attributes
        NcVar timespan_var(ncio.nc->getVar(vname));
        timespan_var.putAtt("units", time_unit.to_cf());
        timespan_var.putAtt("calendar", "365_day");
        timespan_var.putAtt("axis", "T");
        timespan_var.putAtt("long_name", "[start, now] of this time segment");
    }
}
// ----------------------------------------------------------

}
