#ifndef SPSPARSE_ZVECTOR_HPP
#define SPSPARSE_ZVECTOR_HPP

#include <vector>
#include <boost/interprocess/streams/vectorstream.hpp>
#include <boost/endian/conversion.hpp>
#include <boost/enum.hpp>
#include <zstr.hpp>        // https://github.com/mateidavid/zstr
#include <prettyprint.hpp>

namespace spsparse {

// ==============================================================

/** vector-like class, but wraps a std::vector, allows multiple
    boost::interprocess:basic_ivectorstream to read from vector at
    once.

NOTE: This is to the C++11 standard */

template<
    class T,
    class Allocator = std::allocator<T>
>
class _pvector {
    std::vector<T,Allocator> dummy;
    std::vector<T,Allocator> *vec;

    typedef std::vector<T,Allocator> wrapped;
public:
    // Member types (type traits)
    typedef typename wrapped::value_type value_type;
    typedef typename wrapped::allocator_type allocator_type;
    typedef typename wrapped::size_type size_type;
    typedef typename wrapped::difference_type difference_type;
    typedef typename wrapped::reference reference;
    typedef typename wrapped::const_reference const_reference;
    typedef typename wrapped::pointer pointer;
    typedef typename wrapped::const_pointer const_pointer;
    typedef typename wrapped::iterator iterator;
    typedef typename wrapped::const_iterator const_iterator;
    typedef typename wrapped::reverse_iterator reverse_iterator;
    typedef typename wrapped::const_reverse_iterator const_reverse_iterator;

    // Member functions (only those we need)
    _pvector() : vec(&dummy) {}
    _pvector(std::vector<T,Allocator> *_vec) : vec(_vec) {}

    void resize( size_type count )
        { vec->resize(count); }
    void swap( _pvector<T,Allocator>& other )
        { std::swap(vec, other.vec); }
    void reserve( size_type new_cap )
        { vec->reserve(new_cap); }
    void clear()
        { vec->clear(); }
    bool empty() const
        { return vec->empty(); }
    reference operator[]( size_type pos )
        { return vec->operator[](pos); }
    const_reference operator[]( size_type pos ) const
        { return vec->operator[](pos); }
    size_type size() const
        { return vec->size(); }
    size_type capacity() const
        { return vec->capacity(); }
    void push_back( const T& value )
        {vec->push_back(value); }
    void push_back( T&& value )
        { vec->push_back(std::move(value)); }
};

// ==============================================================



enum class ZVAlgo {PLAIN, DIFFS};

// ---------------------------------------------------------
// Yields an integer type, of the same size as the template type

template<class TypeT>
struct ToIntType {
    typedef TypeT int_type;
};

template<>
struct ToIntType<float>
{
    typedef uint32_t int_type;
};

template<>
struct ToIntType<double>
{
    typedef uint64_t int_type;
};
// ---------------------------------------------------------

namespace vaccum {
template<class ValueT, int RANK>
class _ZVector {
public:
    typedef ValueT val_type;

private:
    typedef typename ToIntType<ValueT>::int_type int_type;

    std::vector<char> &zbuf;
    ZVAlgo algo;
    std::array<ValueT,RANK> last_raws;

    boost::interprocess::basic_ovectorstream<std::vector<char>> os;
    zstr::ostream zos;

public:
    _ZVector(std::vector<char> &_zbuf, ZVAlgo _algo);
    void add(std::array<ValueT,RANK> const &raws);
    ~_ZVector();
};

/** Wrap the guts in a std::unique_ptr<> so ZVector is fully C++11
    compliant (and moveable) */
template<class ValueT, int RANK>
class ZVector {
    std::unique_ptr<_ZVector<ValueT,RANK>> self;
public:
    typedef ValueT val_type;

    ZVector(std::vector<char> &_zbuf, ZVAlgo _algo) :
        self(new _ZVector<ValueT,RANK>(_zbuf, _algo)) {}
    void add(std::array<ValueT,RANK> const &raws)
        { self->add(raws); }
};
// -------------------------------------------------------------


/** Compressed version of std::vector */
template<class ValueT, int RANK>
_ZVector<ValueT,RANK>::
    _ZVector(std::vector<char> &_zbuf, ZVAlgo _algo)
        : zbuf(_zbuf),
        algo(_algo),
        os(std::ios_base::out | std::ios_base::binary),
        zos(os, zstr::ostreambuf::default_buff_size, Z_DEFAULT_COMPRESSION)
    {
        // Write the algo to the stream!
        int ialgo = boost::endian::native_to_big((int)algo);
        zos.write((char *)&ialgo, sizeof(ialgo));

        for (int i=0; i<RANK; ++i) last_raws[i] = 0;
    }

template<class ValueT, int RANK>
void _ZVector<ValueT,RANK>::
    add(std::array<ValueT,RANK> const &raws)
    {

        std::array<ValueT,RANK> vals;
        int_type * const ivals((int_type *)&vals[0]);    // Alternative view into vals

        switch(algo) {
            case ZVAlgo::PLAIN :
                for (int i=0; i<RANK; ++i) {
                    vals[i] = raws[i];
                    boost::endian::native_to_big_inplace(ivals[i]);
                }
            break;
            case ZVAlgo::DIFFS :
                for (int i=0; i<RANK; ++i) {
                    vals[i] = raws[i] - last_raws[i];
                    boost::endian::native_to_big_inplace(ivals[i]);
                    last_raws[i] = raws[i];
                }
            break;
        }
        zos.write((char *)ivals, sizeof(int_type)*RANK);
    }

template<class ValueT, int RANK>
_ZVector<ValueT,RANK>::
    ~_ZVector()
    {
        zos.flush();
        // Make output available to the user
        os.swap_vector(zbuf);
    }

// -------------------------------------------------------------

}    // namespace vaccum
// ============================================================================


// --------------------------------------------------------------------------
namespace vgen {

template<class ValueT, int RANK>
class _ZVector {
    typedef ValueT val_type;

private:
    typedef typename ToIntType<ValueT>::int_type int_type;

    _pvector<char> pzbuf;

    std::array<ValueT,RANK> cur_raws;
    ZVAlgo algo;

    boost::interprocess::basic_ivectorstream<_pvector<char>> is;
    zstr::istream zis;

public:
    explicit _ZVector(std::vector<char> &_zbuf);

    bool operator++();

    std::array<ValueT,RANK> const &operator*() const
    {
        return cur_raws;
    }

    ~_ZVector();

};

// -------------------------------------------------------------
/** Wrap the guts in a std::unique_ptr<> so ZVector is fully C++11
    compliant (and moveable) */
template<class ValueT, int RANK>
class ZVector {
    std::unique_ptr<_ZVector<ValueT,RANK>> self;
public:
    typedef ValueT val_type;

    explicit ZVector(std::vector<char> &_zbuf) :
        self(new _ZVector<ValueT,RANK>(_zbuf)) {}

    bool operator++()
        { return self->operator++(); }

    std::array<ValueT,RANK> const &operator*() const
        { return self->operator*(); }
};
// -------------------------------------------------------------

template<class ValueT, int RANK>
_ZVector<ValueT,RANK>::
    _ZVector(std::vector<char> &_zbuf)    // GENERATOR
        : pzbuf(&_zbuf),
        is(std::ios_base::out | std::ios_base::binary),
        zis(is)
    {
        is.swap_vector(pzbuf);

        // Read the algo from the stream
        int ialgo;
        zis.read((char *)&ialgo, sizeof(ialgo));
        boost::endian::big_to_native_inplace(ialgo);
        algo = (ZVAlgo)ialgo;

        for (int i=0; i<RANK; ++i) cur_raws[i] = 0;
    }


template<class ValueT, int RANK>
bool _ZVector<ValueT,RANK>::
    operator++()
    {
        std::array<ValueT,RANK> vals;
        int_type * const ivals((int_type *)&vals[0]);    // Alternative view into vals

        // Read a value
        zis.read((char *)ivals, sizeof(int_type)*RANK);
        if (zis.eof()) return false;

        // Do difference encoding
        switch(algo) {
            case ZVAlgo::PLAIN :
                for (int i=0; i<RANK; ++i) {
                    boost::endian::big_to_native_inplace(ivals[i]);
                    cur_raws[i] = vals[i];
                }
            break;
            case ZVAlgo::DIFFS :
                for (int i=0; i<RANK; ++i) {
                    boost::endian::big_to_native_inplace(ivals[i]);
                    cur_raws[i] += vals[i];
                }
            break;
        }
        return true;
    }

template<class ValueT, int RANK>
_ZVector<ValueT,RANK>::
    ~_ZVector()
    {
        // Get our original vector back!
        // (We don't need this with _pvector
        // // is.swap_vector(zbuf);
    }

}    // namespace spsparse::vgen


}    // Namespace spsparse
#endif    // guard
