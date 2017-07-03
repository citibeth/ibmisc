#ifndef IBMISC_FORTRANIO_HPP
#define IBMISC_FORTRANIO_HPP

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <blitz/array.h>
#include <ibmisc/error.hpp>
#include <ibmisc/endian.hpp>

namespace ibmisc {
namespace fortran {

struct UnformattedInput {
    std::ifstream fin;
    Endian const endian;

    UnformattedInput(std::string const &fname, ibmisc::Endian _endian) :
        fin(fname, std::ios::binary | std::ios::in),
        endian(_endian) {}

};


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
    size_t const len;    // Length to read, in bytes

    BufSpec(size_t const _len) : len(_len) {}
    virtual ~BufSpec() {}
    virtual void read(UnformattedInput &infile) = 0;
};

struct SimpleBufSpec : public BufSpec {
    char * const buf;
    int item_size;    // 1, 2, 4, 8
    long nitem;

    SimpleBufSpec(char *_buf, int item_size, long nitem);
    void read(UnformattedInput &infile);
};

class EndR {};

extern EndR endr;

/** Wrap double arrays in this when we want to read into a float. */
template<class SrcT, class DestT, int RANK>
class BlitzCast : public BufSpec
{
    blitz::Array<DestT, RANK> dest;

public:
    BlitzCast(blitz::Array<DestT, RANK> &_dest);

    void read(UnformattedInput &infile);
};


template<class SrcT, class DestT, int RANK>
BlitzCast<SrcT, DestT, RANK>::BlitzCast(blitz::Array<DestT, RANK> &_dest) :
    BufSpec(_dest.size() * sizeof(SrcT)),
    dest(_dest)
{
    if (!dest.isStorageContiguous()) (*ibmisc_error)(-1,
        "Storage must be contiguous");
}


template<class SrcT, class DestT, int RANK>
void BlitzCast<SrcT, DestT, RANK>::read(UnformattedInput &infile)
{
    blitz::Array<SrcT, 1> src(dest.size());
    char * const buf = (char *)src.data();
    infile.fin.read(buf, len);

     // Swap for endian
    endian_to_native(buf, sizeof(SrcT), dest.size(), infile.endian);

    auto srci(src.begin());
    auto desti(dest.begin());
    for (; srci != src.end(); ++srci, ++desti) {
        *desti = (DestT)*srci;
    }
}


template<class SrcT, class DestT, int RANK>
std::unique_ptr<BufSpec> blitz_cast(blitz::Array<DestT, RANK> &dest)
    { return std::unique_ptr<BufSpec>(new BlitzCast<SrcT,DestT,RANK>(dest)); }

// ---------------------------------------------------------------
class read {
    UnformattedInput *infile;
    std::vector<std::unique_ptr<BufSpec>> specs;

    read &add_simple(char *buf, int item_size, long nitem)
    {
        specs.push_back(std::unique_ptr<BufSpec>(new SimpleBufSpec(buf, item_size, nitem)));
        return *this;
    }


public:
    read(UnformattedInput &_infile) : infile(&_infile) {}

    read &operator>>(std::unique_ptr<BufSpec> &&specp)
    {
        specs.push_back(std::move(specp));
        return *this;
    }

    template<class TypeT, int RANK>
    read &operator>>(blitz::Array<TypeT,RANK> &arr)
    {
        return add_simple((char *)arr.data(), sizeof(TypeT), arr.size());
    }


    template<class TypeT, size_t SIZE>
    read &operator>>(std::array<TypeT, SIZE> &arr)
    {
        return add_simple((char *)&arr[0], sizeof(TypeT), arr.size());
    }

    read &operator>>(float &val)
    {
        return add_simple((char *)&val, sizeof(float), 1);
    }

    read &operator>>(int &val)
    {
        return add_simple((char *)&val, sizeof(int), 1);
    }

    read &operator>>(double &val)
    {
        return add_simple((char *)&val, sizeof(double), 1);
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
