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

#include <algorithm>
#include <string>
#include <vector>
#include <cstdarg>
#include <iostream>
#include <sstream>

namespace ibmisc {

std::string string_printf(const std::string& format, ...);

inline void toupper(std::string &str) {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

inline void tolower(std::string &str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

inline bool ends_with(std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

inline bool starts_with(std::string const &fullString, std::string const &prefix)
{
    if (fullString.length() >= prefix.length()) {
        return (0 == fullString.compare(0, prefix.length(), prefix));
    } else {
        return false;
    }
}

bool replace(std::string& str, const std::string& from, const std::string& to);

// Convert Fortran to C++ string
inline std::string f_to_cpp(char *fstr, size_t len)
{
    // Starting from last character, find first non-blank
    char *c;
    for (c=fstr+len-1; *c == ' ' && c > fstr; --c) ;

    // Copy to a C++ string
    return std::string (fstr, c+1-fstr);
}


/** Wrap a paragraph of text, splitting on word boundaries.
https://www.rosettacode.org/wiki/Word_wrap#C.2B.2B */
extern std::string wrap(std::string const &text, size_t line_length = 72);

}
