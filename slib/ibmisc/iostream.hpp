#ifndef IBMISC_IOSTREAM
#define IBMISC_IOSTREAM

#include <iostream>

/** Write std::array to ostream. */
template <class T, std::size_t N>
ostream& operator<<(ostream& o, const std::array<T, N>& arr)
{
    o << "[";
    copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, " "));
    0 << "]";
    return o;
}

#endif
