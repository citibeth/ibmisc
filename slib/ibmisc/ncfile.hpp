#ifndef IBMISC_NCFILE_HPP
#define IBMISC_NCFILE_HPP

#include <ibmisc/netcdf.hpp>
#include <ibmisc/datetime.hpp>

namespace ibmisc {


/** For files with a single timepoint: Defines start and end timepoints.
@param timespan Start and end of the current timestep's time.
@param time_units Units of time in NetCDF format (eg: "seconds since 2001-01-01") */
extern void ncio_timespan(NcIO &ncio, std::array<double,2> &timespan, ibmisc::TimeUnit const &time_unit, std::string const &vname);


}
#endif
