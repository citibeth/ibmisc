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

#include <ibmisc/string.hpp>

namespace ibmisc {

/** Like sprintf(), but returns output in a std::string.  Uses same
parameters as printf().
@see http://codereview.stackexchange.com/questions/52522/mimic-sprintf-with-stdstring-output */
std::string string_printf(const std::string& format, ...)
{
    static const int initial_buf_size = 100;
    va_list arglist;
    va_start(arglist, format);
    char buf1[initial_buf_size];
    const int len = vsnprintf(buf1,initial_buf_size,format.c_str(), arglist) + 1;
    va_end(arglist);

    if(len<initial_buf_size){
        return buf1;
    } else {
        char buf2[len];
        va_start(arglist,format);
        vsnprintf(buf2,len,format.c_str(),arglist);
        va_end(arglist);
        return buf2;
    }
}


// See: http://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}


}
