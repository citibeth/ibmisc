#include <cstdio>
#include <cstdarg>
#include <exception>
#include <ibmisc/ibmisc.hpp>
#ifdef USE_EVERYTRACE
#include <everytrace.h>
#endif

namespace ibmisc {

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
	throw ibmisc::Exception();
//	exit(-1);
}

error_ptr ibmisc_error = &default_error;

} 	// Namespace
