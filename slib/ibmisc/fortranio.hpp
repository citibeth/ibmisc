#include <string>
#include <vector>
#include <array>

#include <blitz/array.h>
#include <ibmisc/error.hpp>

namespace ibmisc {
namespace fortran {

/** This namespace provides FORTRAN-like functions to read unformatted
    record-based binary Fortran files.  The idea is NOT to be more
    clever than Fortran; but rather, to faithfully replicate the way
    Fortran does it.

Example:

    #include <iostream>
    #include <fstream>

    std::array<char, 80> str0, str1;
    blitz::Array<int, 1> vals(17);

    // Open the file for reading, using standard C++ oeners
    std::ifstream fin("sample_fortranio", ios::binary | ios::in);

    // Read a record into the variables str0 and vals
    fortran::read(fin) >> str0 >> vals >> fortran::endr;

    // Trim (remove trailing whitespace from) the string we read
    str0 = fortran::trim(str0);

    // Skip a record    
    fortran::read(fin) >> fortran::endr;

    // Read a record into the variables str1, vals and str0
    fortran::read(fin) >> str1 >> vals >> str0 >> fortran::endr;
*/

struct BufSpec {
    char * const buf;
    size_t const len;

    BufSpec(char *_buf, size_t _len) : buf(_buf), len(_len) {}
};

class EndR {};

extern EndR endr;

class read {
    std::istream *infile;
    std::vector<BufSpec> specs;
public:
    read(std::istream &_infile) : infile(&_infile) {}

    template<class TypeT, int RANK>
    read &operator>>(blitz::Array<TypeT,RANK> &arr)
    {
        specs.push_back(BufSpec((char *)arr.data(), sizeof(TypeT) * arr.size()));
        return *this;
    }

    template<class TypeT, size_t SIZE>
    read &operator>>(std::array<TypeT, SIZE> &arr)
    {
        specs.push_back(BufSpec((char *)&arr[0], sizeof(TypeT) * arr.size()));
        return *this;
    }

    read &operator>>(float &val)
    {
        specs.push_back(BufSpec((char *)&val, sizeof(float)));
        return *this;
    }

    read &operator>>(int &val)
    {
        specs.push_back(BufSpec((char *)&val, sizeof(int)));
        return *this;
    }

    read &operator>>(double &val)
    {
        specs.push_back(BufSpec((char *)&val, sizeof(double)));
        return *this;
    }



    // Do the read!
    void operator>>(EndR const &endr);

};

// -------------------------------------------------------
/** Trims trailing whitespace.  Returns a new string, because
that is how Fortran's TRIM() works. */
template<size_t SIZE>
std::string trim(std::array<char, SIZE> const &fstr);

template<size_t SIZE>
std::string trim(std::array<char, SIZE> const &fstr)
{
    size_t i=SIZE-1;
    for (; i >= 0; --i) if (fstr[i] != ' ') break;
    return std::string(&fstr[0], i+1);
}


}}    // ibmisc::fortran
