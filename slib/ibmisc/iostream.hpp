#ifndef IBMISC_IOSTREAM
#define IBMISC_IOSTREAM

#include <iostream>

/** Write std::array to ostream. */
template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr)
{
    o << "[";
    copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
    o << "]";
    return o;
}

/** Write std::vector to ostream. */
template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& arr)
{
    o << "[";
    copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
    o << "]";
    return o;
}

#endif
