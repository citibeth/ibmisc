#ifndef IBMISC_NCFILE_HPP
#define IBMISC_NCFILE_HPP

#include <ibmisc/netcdf.hpp>

namespace ibmisc {


#if 0
/** For files with multiple timepoints: Defines a time0 variable
    showing start time, and a time variable for each time point in the
    file.
@param time0 Time of the start of simulation
@param time_units Units of time in NetCDF format (eg: "seconds since 2001-01-01") */
extern void ncdefine_time0(NcIO &ncio, double time0, std::string const &time_units, std::string const &vname_base);
#endif

/** For files with a single timepoint: Defines start and end timepoints.
@param timespan Start and end of the current timestep's time.
@param time_units Units of time in NetCDF format (eg: "seconds since 2001-01-01") */
extern void ncio_timespan(NcIO &ncio, std::array<double,2> &timespan, std::string const &time_units, std::string const &vname_base);


}
#endif
