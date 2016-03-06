#include <cstdio>
#include <cstdarg>
#include <exception>
#include <stdexcept>
#include <ibmisc/ibmisc.hpp>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

namespace ibmisc {


// http://codereview.stackexchange.com/questions/52522/mimic-sprintf-with-stdstring-output

std::string vsprintf(const char* format, std::va_list args)
{
    va_list tmp_args; //unfortunately you cannot consume a va_list twice
    va_copy(tmp_args, args); //so we have to copy it
    const int required_len = vsnprintf(nullptr, 0, format, tmp_args) + 1;
    va_end(tmp_args);

    std::string buf(required_len, '\0');
    if (std::vsnprintf(&buf[0], buf.size(), format, args) < 0) {
        (*ibmisc_error)(-1, "string_vsprintf encoding error");
    }
    return buf;
}

std::string sprintf(const char* format, ...)
{
    std::va_list args;
    va_start(args, format);
    std::string str{ibmisc::vsprintf(format, args)};
    va_end(args);
    return str;
}
// --------------------------------------------------------


void exception_error(int retcode, const char *format, ...)
{
    std::va_list args;
    va_start(args, format);
    std::string str{ibmisc::vsprintf(format, args)};
    va_end(args);
    throw std::runtime_error(std::move(str));
}

#ifdef USE_EVERYTRACE
void everytrace_error(int retcode, const char *format, ...)
{
    va_list arglist;

    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);
    fprintf(stderr, "\n");

    everytrace_exit(retcode);
}
error_ptr ibmisc_error = &everytrace_error;
#else
error_ptr ibmisc_error = &exception_error;

#endif



}   // Namespace
