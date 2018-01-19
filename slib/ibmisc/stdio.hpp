#ifndef IBMISC_STDIO_HPP
#define IBMISC_STDIO_HPP

#include <string>
#include <cstdarg>

namespace ibmisc {

// http://codereview.stackexchange.com/questions/52522/mimic-sprintf-with-stdstring-output

extern std::string vstrprintf(const char* format, std::va_list args);

extern std::string strprintf(const char* format, ...);

}    // namespace

#endif
