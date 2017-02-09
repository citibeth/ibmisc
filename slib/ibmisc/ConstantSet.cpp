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

#include <cstring>
#include <ibmisc/ConstantSet.hpp>
#include <ibmisc/string.hpp>

using namespace netCDF;

namespace ibmisc {

int ConstantSet::add(
    std::string const &name, std::string const &units,
    std::string const &description)
{
    index.insert(name);
    data.push_back(Data(name, units, description));
    return index.size()-1;
}

int ConstantSet::set(
    std::string const &name,
    double val,
    std::string const &units,
    std::string const &description)
{
    int ix = add(name, units, description);
    data[ix].val = val;
    return ix;
}

/** Copy a constant's value and units to another constant
@return Index of the new constant. */
int ConstantSet::copy(std::string const &dst_name,
    std::string const &src_name,
    std::string const &description)
{
    int src_ix = index.at(src_name);
    int dst_ix = add(dst_name, data[src_ix].units, description);
    data[dst_ix].val = data[src_ix].val;

    return dst_ix;
}


double ConstantSet::get_as(std::string const &name,
    UTUnit const &units) const
{
    int src_ix = index.at(name);
    UTUnit usrc(ut_system->parse(data[src_ix].units));

    try {
        CVConverter cv(usrc, units);
        double ret = cv.convert(data[src_ix].val);
        double src_val = (*this)[src_ix].val;
        if (ret != src_val) {
            printf("ConstantSet: Converting %s: %g [%s] --> %g [%s]\n", name.c_str(), src_val, usrc.c_str(), ret, units.c_str());
        }
        return ret;
    } catch(const std::exception &ex) {
        (*ibmisc_error)(-1,
            "Exception in ConstantSet::get_as(%s, %s)\n", name.c_str(), units.c_str());
    }
    
}

double ConstantSet::get_as(std::string const &name,
    std::string const &sunits) const
{
    UTUnit units(ut_system->parse(sunits));
    return get_as(name, units);
}

// =======================================================
void ConstantSet::ncio(NcIO &ncio, std::string const &vname)
{
    auto constants_v = get_or_add_var(ncio, vname, "int64", {});

    // Get a list of the names of our constants
    if (ncio.rw == 'w') {
        // Store the constants as attributes
        for (size_t i=0; i<size(); ++i) {
            constants_v.putAtt(data[i].name, ncDouble, data[i].val);
            constants_v.putAtt(data[i].name + "_units", data[i].units);
            constants_v.putAtt(data[i].name + "_description", data[i].description);
        }
    } else {
        auto atts(constants_v.getAtts());   // Get all atts for this var
        for (auto ii=atts.begin(); ii != atts.end(); ++ii) {

            // Loop through the MAIN attribute for each constant
            std::string const &name(ii->first);
            NcAtt &att(ii->second);
            if (att.getType() != ncDouble) continue;

            // Get its value
            if (att.getAttLength() != 1) (*ibmisc_error)(-1,
                "Attribute %s must have length 1", name.c_str());
            double val;
            att.getValues(&val);

            // Look up other attributes
            std::string units, description;
            get_or_put_att(constants_v, 'r', name + "_units", units);
            get_or_put_att(constants_v, 'r', name + "_description", description);

            // Create the constant
            set(name, val, units, description);
        }
    }
}

}
