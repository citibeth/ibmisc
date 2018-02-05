#ifndef IBMISC_FORTRANIO_HPP
#define IBMISC_FORTRANIO_HPP

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <ibmisc/blitz.hpp>
#include <ibmisc/error.hpp>
#include <ibmisc/endian.hpp>
#include <ibmisc/memory.hpp>

namespace ibmisc {
namespace fortran {

struct UnformattedInput {
    std::ifstream fin;
    Endian const endian;
    TmpAlloc tmp;

    UnformattedInput(std::string const &fname, ibmisc::Endian _endian) :
        fin(fname, std::ios::binary | std::ios::in),
        endian(_endian) {}

    void close() { fin.close(); }

    bool eof() const { return fin.eof(); }

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
    long nbytes;    // Length to read, in bytes (-1 = wildcard, read all available)
    bool const wildcard;    // Was this constructed with len<0 (wildcard matching length)?

    BufSpec(size_t const _nbytes) : nbytes(_nbytes), wildcard(nbytes<0) {}

    virtual ~BufSpec() {}
    // Call on wildcards
    virtual void set_nbytes(long _nbytes)
        { (*ibmisc_error)(-1, "BufSpec does not support set_nbytes()"); }
    virtual void read(UnformattedInput &infile) = 0;
};
class EndR {};
extern EndR endr;

struct SimpleBufSpec : public BufSpec {
    char * const buf;
    int item_size;    // 1, 2, 4, 8
    long nitem;

    SimpleBufSpec(char *_buf, int item_size, long nitem);
    void read(UnformattedInput &infile);
};

/** This class is here to create a bunch of implicit conversions */
struct BufSpecPtr : public std::unique_ptr<BufSpec>
{
    typedef std::unique_ptr<BufSpec> super;

    BufSpecPtr(BufSpec *spec) : super(spec) {}

    template<class TypeT, int RANK>
    BufSpecPtr(blitz::Array<TypeT,RANK> &buf);

    template<class TypeT, size_t SIZE>
    BufSpecPtr(std::array<TypeT, SIZE> &arr)
        : super(new SimpleBufSpec((char *)&arr[0], sizeof(TypeT), arr.size())) {}

    BufSpecPtr(float &val)
        : super(new SimpleBufSpec((char *)&val, sizeof(float), 1)) {}
    BufSpecPtr(int &val)
        : super(new SimpleBufSpec((char *)&val, sizeof(int), 1)) {}
    BufSpecPtr(double &val)
        : super(new SimpleBufSpec((char *)&val, sizeof(double), 1)) {}

};
// -----------------------------------------------------------------------------
template<int RANK>
struct Shape {
    std::array<std::string, RANK> sshape;
    std::array<int,RANK> shape;

    Shape(
        std::array<std::string, RANK> &&_sshape,
        std::array<int,RANK> &&_shape)
        : sshape(std::move(_sshape)), shape(std::move(_shape))
    {}

    long size() const {
        long ret = 1;
        for (int i=0; i<RANK; ++i) ret *= shape[i];    
        return ret;
    }
};

// ---------------------------------------------------------------------
/** Reads 1-D arrays */
template<class TypeT, int RANK>
struct ArrayBuf : public BufSpec {
    blitz::Array<TypeT,RANK> *buf;
    Shape<RANK> const **buf_shape;
    blitz::vector<Shape<RANK>> const *stdshapes;        // Possible shapes we might infer from nbytes
    blitz::GeneralArrayStorage<RANK> const storage;

    /** User-supplied pre-allocated buffer */
    ArrayBuf(
        blitz::Array<TypeT,RANK> *_buf)
    : BufSpec(_buf->size() * sizeof(TypeT)), buf(_buf) {}

    /** Wildcard: will allocate the array */
    ArrayBuf(
        blitz::Array<TypeT,RANK> *_buf,
        Shape<RANK> const **_buf_shape,
        blitz::vector<Shape<RANK>> const *_stdshapes,
        blitz::GeneralArrayStorage<RANK> const &_storage = blitz::fortranArray)
    : BufSpec(-1), buf(_buf), buf_shape(_buf_shape), stdshapes(_stdshapes), storage(_storage)
    {}

    /** Wildcard: will allocate the array, no user-supplied buffer at all. */
    ArrayBuf(
        blitz::vector<Shape<RANK>> const *_stdshapes,
        blitz::GeneralArrayStorage<RANK> const &_storage = blitz::fortranArray)
    : BufSpec(-1), buf(0), stdshapes(_stdshapes), storage(_storage)
    {}

    void set_nbytes(long _nbytes)
    {
        nbytes = _nbytes;
    }

    virtual void read(UnformattedInput &infile);
};

template<class TypeT, int RANK>
void ArrayBuf<TypeT,RANK>::read(UnformattedInput &infile)
{
    // Allocate the array, if this is a wildcard reader...
    if (wildcard) {
        blitz::TinyVector<int,RANK> shape;    // Dimensions we should allocate
        if (RANK == 1) {
            shape[0] = nbytes / sizeof(TypeT);
        } else {
            for (auto sh(stdshapes->begin()); ; ++sh) {
                if (sh == stdshapes->end()) (*ibmisc_error)(-1,
                    "Cannot find std shape for nbytes=%ld (sizeof=%ld)", nbytes, sizeof(TypeT));

                if (sh->size() * sizeof(TypeT) == nbytes) {
                    for (int i=0; i<RANK; ++i) shape[i] = sh->shape[i];
                    if (buf_shape) *buf_shape = &*sh;
                    break;
                }
            }
        }
        if (!buf) {
            buf = infile.tmp.newptr<blitz::Array<TypeT,RANK>>();
        }
        buf->reference(blitz::Array<TypeT,RANK>(shape, storage));
    }

    // Read it!
    SimpleBufSpec sbuf((char *)buf->data(), sizeof(TypeT), buf->size());
    sbuf.read(infile);
}
        

/** User-supplied pre-allocated buffer */
template<class TypeT, int RANK>
BufSpecPtr::BufSpecPtr(blitz::Array<TypeT,RANK> &buf)
    : super(new ArrayBuf<TypeT,RANK>(&buf))
    {}

/** Wildcard: will allocate the array */
template<class TypeT, int RANK>
BufSpecPtr allocate(
    blitz::Array<TypeT,RANK> &_buf,
    Shape<RANK> const *&_buf_shape,
    blitz::vector<Shape<RANK>> const &_stdshapes,
    blitz::GeneralArrayStorage<RANK> const &_storage = blitz::fortranArray)
{ return BufSpecPtr(new ArrayBuf<TypeT,RANK>(&_buf, &_buf_shape, &_stdshapes, _storage)); }

template<class TypeT, int RANK>
BufSpecPtr allocate(
    blitz::Array<TypeT,RANK> &_buf,
    blitz::vector<Shape<RANK>> const &_stdshapes,
    blitz::GeneralArrayStorage<RANK> const &_storage = blitz::fortranArray)
{ return BufSpecPtr(new ArrayBuf<TypeT,RANK>(&_buf, NULL, &_stdshapes, _storage)); }

// ============================================================================
// ---------------------------------------------------------------------
/** Wrap double arrays in this when we want to read into a float. */
template<class SrcT, class DestT, int RANK>
class BlitzCast : public BufSpec
{
    TmpAlloc tmp;
    blitz::Array<DestT, RANK> *dest;
    std::unique_ptr<ArrayBuf<SrcT,RANK>> subbuf;


public:
    /** Cast into our pre-allocated array */
    BlitzCast(
        blitz::Array<DestT,RANK> *_dest)
    : BufSpec(sizeof(SrcT) * _dest->size()),
        dest(_dest),
        subbuf(new ArrayBuf<SrcT,RANK>(tmp.newptr<blitz::Array<SrcT,RANK>>(_dest->shape())))
    {}

    /** Cast into a wildcard */
    BlitzCast(
        blitz::Array<DestT,RANK> *_dest,
        blitz::vector<Shape<RANK>> const *_stdshapes,
        blitz::GeneralArrayStorage<RANK> const &_storage = blitz::fortranArray)
    : BufSpec(-1), dest(_dest),
        subbuf(new ArrayBuf<SrcT,RANK>(_stdshapes, _storage)) {}

    void set_nbytes(long _nbytes)
        { subbuf->set_nbytes(_nbytes); }

    void read(UnformattedInput &infile);
};


template<class SrcT, class DestT, int RANK>
void BlitzCast<SrcT, DestT, RANK>::read(UnformattedInput &infile)
{
    subbuf->read(infile);

    if (subbuf->wildcard) {
        dest->reference(blitz::Array<DestT,RANK>(subbuf->buf->shape(), subbuf->storage));
    } else {
        check_dimensions("array", *dest, tiny_to_vector(subbuf->buf->shape()));
    }

    // Copy the data over
    auto srci(subbuf->buf->begin());
    auto desti(dest->begin());
    for (; srci != subbuf->buf->end(); ++srci, ++desti) {
        *desti = (DestT)*srci;
    }
}




/** Cast into our pre-allocated array */
template<class SrcT, class DestT, int RANK>
BufSpecPtr blitz_cast(blitz::Array<DestT,RANK> &_dest)
    { return BufSpecPtr(new BlitzCast<SrcT,DestT,RANK>(&_dest)); }

/** Cast into a wildcard */
template<class SrcT, class DestT, int RANK>
BufSpecPtr blitz_cast(
    blitz::Array<DestT,RANK> &_dest,
    blitz::vector<Shape<RANK>> const *_stdshapes,
    blitz::GeneralArrayStorage<RANK> const &_storage = blitz::fortranArray)
    { return BufSpecPtr(new BlitzCast<SrcT,DestT,RANK>(&_dest, _stdshapes, _storage)); }

// ---------------------------------------------------------------
class read {
    UnformattedInput *infile;
    std::vector<std::unique_ptr<BufSpec>> specs;

public:
    read(UnformattedInput &_infile) : infile(&_infile) {}

    read &operator>>(BufSpecPtr &&specp)
    {
        specs.push_back(std::move(specp));
        return *this;
    }


    // Do the read! (TODO: Do this with a destructor)
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
