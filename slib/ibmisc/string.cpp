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

}