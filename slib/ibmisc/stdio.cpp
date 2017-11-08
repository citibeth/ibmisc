#include <ibmisc/stdio.hpp>
#include <ibmisc/error.hpp>

namespace ibmisc {

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

}    // namespace
