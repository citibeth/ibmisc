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

#pragma once

#include <string>
#include <boost/enum.hpp>
#include <ibmisc/ibmisc.hpp>

namespace ibmisc {

template<class T>
inline T parse_enum(std::string const &str) {
    auto ret = T::get_by_name(str.c_str());
    if (!ret) {
        (*ibmisc_error)(-1,
            "Error converting from string '%s' for boost::enum type %s\n", str.c_str(), typeid(T).name());
    }
    return *ret;
}

}
