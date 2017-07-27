#ifndef IBMISC_ENDIAN_HPP
#define IBMISC_ENDIAN_HPP

#include <iostream>
#include <boost/lexical_cast.hpp>

namespace ibmisc {

enum class Endian {LITTLE, BIG};

/** Convert an entire array of big-endian items to native-endian format */
void big_to_native(char *buf, int const item_size, long nitem);

/** Convert an entire array of big-endian items to native-endian format */
void little_to_native(char *buf, int const item_size, long nitem);

/** Convert big or little to native */
void endian_to_native(char *buf, int const item_size, long nitem, Endian endian);

} // namespace

// =================================================================
// Used for boost::program_options

namespace boost {

template<>
ibmisc::Endian lexical_cast<ibmisc::Endian, std::string>(std::string const &token);

template<>
std::string lexical_cast<std::string, ibmisc::Endian>(ibmisc::Endian const &endian);

}    // namespace boost

#endif    // guard
