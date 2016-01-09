/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <ibmisc/Proj2.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/ibmisc.hpp>

namespace ibmisc {

Proj2::Proj2(std::string const &_sproj, Direction _direction) :
sproj(_sproj), direction(_direction)
	{ realize(); }


Proj2::Proj2(Proj2 const &rhs) :
sproj(rhs.sproj), direction(rhs.direction)
	{ realize(); }


void Proj2::realize()
{
	if (sproj == "") return;

	_proj = Proj(sproj);
	_llproj = _proj.latlong_from_proj();
}

/** Transforms a single coordinate pair
@param src Source coordinate system
@param dest Destination coordinate system.
@param x0 Source x (or longitude) coordinate (radians)
@param y0 Source y (or latitude) coordinate (radians)
@param x1 Destination x (or longitude) coordinate (radians)
@param y1 Destination y (or latitude) coordinate (radians) */
int Proj2::transform(double x0, double y0, double &x1, double &y1) const
{
	if (!is_valid()) {
		x1 = x0;
		y1 = y0;
		return 0;
	}

	if (direction == Direction::XY2LL) {
		int ret = ibmisc::transform(_proj, _llproj, x0, y0, x1, y1);
		x1 *= R2D;
		y1 *= R2D;
		return ret;
	}

	x0 *= D2R;
	y0 *= D2R;
	return ibmisc::transform(_llproj, _proj, x0, y0, x1, y1);
}


#ifdef USE_NETCDF

void ncio_proj2(
	ibmisc::NcIO &ncio,
	std::string const &vname,
	Proj2 &proj,
	std::string const &attrname)
{
	netCDF::NcVar ncvar = ncio.nc->getVar(vname);
	if (ncio.rw == 'w') {
		// --------- WRITE
		if (proj.is_valid()) {
			ncvar.putAtt(attrname, proj.sproj);
			std::string sdir = (proj.direction == Proj2::Direction::XY2LL ? "xy2ll" : "ll2xy");
			ncvar.putAtt((attrname + ".direction"), sdir);
		} else {
			ncvar.putAtt(attrname, "");
			ncvar.putAtt(attrname + ".direction", "");
		}
	} else {
		// ----------- READ
		std::string sproj;
		ncvar.getAtt(attrname).getValues(sproj);
		if (sproj == "") {
			proj.clear();
		} else {
			std::string sdir;
			ncvar.getAtt(attrname + ".direction").getValues(sproj);
			if (sdir == "xy2ll") proj.direction = Proj2::Direction::XY2LL;
			else if (sdir == "ll2xy") proj.direction = Proj2::Direction::LL2XY;
			else (*ibmisc_error)(-1,
				"Unsupported direction \"%s\" from NetCDF file (must be \"xy2ll\" or \"ll2xy\"", sdir.c_str());
		}
	}
}

#endif

}	// namespace ibmisc
