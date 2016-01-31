#ifndef IBMISC_IBMISC_HPP
#define IBMISC_IBMISC_HPP

#include <string>
#include <cstdarg>
#include <stdexcept>

/** @defgroup ibmisc ibmisc.hpp
@brief Basic stuff common to all ibmisc */
namespace ibmisc {

typedef std::runtime_error Exception;

std::string vsprintf(const char* format, std::va_list args);

std::string sprintf(const char* format, ...) __attribute__ ((format (printf, 1, 2)));



// Use this instead.
// http://www.thecodingforums.com/threads/function-pointers-to-printf.317925/
/** @brief Printf-like signature of error handle functions to be used by SpSparse. */
typedef void (*error_ptr) (int retcode, char const *str, ...);

/** @brief Error handler used by IBMisc.  May be changed by user's
main program, to fit into some larger error handling system (eg:
Everytrace).

https://github.com/citibob/everytrace */
extern ibmisc::error_ptr ibmisc_error;

#ifdef USE_EVERYTRACE
extern void everytrace_error(int retcode, const char *format, ...);
#endif
extern void exception_error(int retcode, const char *format, ...);

}	// namespace
/** @} */

#endif // Guard
