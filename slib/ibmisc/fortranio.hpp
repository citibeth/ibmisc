#ifndef IBMISC_FORTRANIO_HPP
#define IBMISC_FORTRANIO_HPP

#include <memory>
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

struct BufSpecI {
    size_t const len;    // Length to read, in bytes

    BufSpecI(size_t const _len) : len(_len) {}
    virtual ~BufSpecI() {}
    virtual void read(std::istream &infile) = 0;
};

struct BufSpec : public BufSpecI {
    char * const buf;

    BufSpec(char *_buf, size_t _len) : BufSpecI(_len), buf(_buf) {}
    void read(std::istream &infile)
    {
        infile.read(buf, len);
    }
};

class EndR {};

extern EndR endr;

/** Wrap double arrays in this when we want to read into a float. */
template<class SrcT, class DestT, int RANK>
class blitz_cast : public BufSpecI
{
    blitz::Array<double, RANK> dest;

public:
    blitz_cast(blitz::Array<DestT, RANK> &dest) :
        BufSpecI(dest.size() * sizeof(SrcT))
    {
        if (!dest.isStorageContiguous()) (*ibmisc_error)(-1,
            "Storage must be contiguous");
    }

    void read(std::istream &infile)
    {
        blitz::Array<float, 1> src(dest.size());

        infile.read((char *)src.data(), len);

        auto srci(src.begin());
        auto desti(dest.begin());
        for (; srci != src.end(); ++srci, ++desti) {
            *desti = (DestT)*srci;
        }
    }
};

class read {
    std::istream *infile;
    std::vector<std::unique_ptr<BufSpec>> specs;

    template<class BufSpecT>
    read &add(BufSpecT &&spec)
    {
        specs.push_back(std::unique_ptr<BufSpec>((BufSpec *)new BufSpecT(std::move(spec))));
        return *this;
    }

public:
    read(std::istream &_infile) : infile(&_infile) {}

    template<class TypeT>
    read &operator>>(TypeT val) {};



    template<class TypeT, int RANK>
    read &operator>>(blitz::Array<TypeT,RANK> &arr)
    {
        return this->add(BufSpec((char *)arr.data(), sizeof(TypeT) * arr.size()));
    }


    template<class TypeT, size_t SIZE>
    read &operator>>(std::array<TypeT, SIZE> &arr)
    {
        return this->add(BufSpec((char *)&arr[0], sizeof(TypeT) * arr.size()));
    }

    read &operator>>(float &val)
    {
        return this->add(BufSpec((char *)&val, sizeof(float)));
    }

    read &operator>>(int &val)
    {
        return this->add(BufSpec((char *)&val, sizeof(int)));
    }

    read &operator>>(double &val)
    {
        return this->add(BufSpec((char *)&val, sizeof(double)));
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

#endif // guard
