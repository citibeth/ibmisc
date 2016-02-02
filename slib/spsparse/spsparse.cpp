#include <functional>
#include <cstdlib>
#include <spsparse/spsparse.hpp>
#include <exception>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif


namespace spsparse {

void default_error(int retcode, const char *format, ...)
{
	va_list arglist;

	va_start(arglist, format);
	vfprintf(stderr, format, arglist);
	va_end(arglist);
	fprintf(stderr, "\n");

#ifdef USE_EVERYTRACE
	everytrace_exit(-1);
#endif
	throw spsparse::Exception();
//	exit(-1);
}

error_ptr spsparse_error = &default_error;

const std::array<int,2> ROW_MAJOR = {0,1};
const std::array<int,2> COL_MAJOR = {1,0};

} 	// Namespace
