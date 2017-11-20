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
void ConstantSet::read_nc(netCDF::NcGroup *nc, std::string const &prefix)
{
    std::multimap<std::string,NcVar> vars(nc->getVars());
    for (auto ii=vars.begin(); ii != vars.end(); ++ii) {
        std::string const &var_name(ii->first);
        if (!starts_with(var_name, prefix)) continue;
        NcVar &var(ii->second);

        std::string const const_name(var_name.substr(prefix.length(), var_name.length() - prefix.length()));


        double value;
        get_or_put_att(var, 'r', "value", &value, 1);
        std::string units;
        get_or_put_att(var, 'r', "units", units);
        std::string description;
        get_or_put_att(var, 'r', "description", description);

        set(const_name, value, units, description);
    }
}

}
