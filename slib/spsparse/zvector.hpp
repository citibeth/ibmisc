#ifndef SPSPARSE_ZVECTOR_HPP
#define SPSPARSE_ZVECTOR_HPP

#include <vector>
#include <boost/interprocess/streams/vectorstream.hpp>
#include <boost/endian/conversion.hpp>
#include <boost/enum.hpp>
#include <zstr.hpp>        // https://github.com/mateidavid/zstr
#include <prettyprint.hpp>

namespace spsparse {

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
class ZVector {
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
    ZVector(std::vector<char> &_zbuf, ZVAlgo _algo);
    void add(std::array<ValueT,RANK> const &raws);
    ~ZVector();
};


template<class ValueT, int RANK>
ZVector<ValueT,RANK>::
    ZVector(std::vector<char> &_zbuf, ZVAlgo _algo)
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
void ZVector<ValueT,RANK>::
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
ZVector<ValueT,RANK>::
    ~ZVector()
    {
        zos.flush();
        // Make output available to the user
        os.swap_vector(zbuf);
    }
}
// ============================================================================



// --------------------------------------------------------------------------
namespace vgen {

template<class ValueT, int RANK>
class ZVector {
    typedef ValueT val_type;

private:
    typedef typename ToIntType<ValueT>::int_type int_type;

    std::vector<char> &zbuf;

    std::array<ValueT,RANK> cur_raws;
    ZVAlgo algo;

    boost::interprocess::basic_ivectorstream<std::vector<char>> is;
    zstr::istream zis;

public:
    explicit ZVector(std::vector<char> &_zbuf);

    bool operator++();

    std::array<ValueT,RANK> const &operator*() const
    {
        return cur_raws;
    }

    ~ZVector();

};


template<class ValueT, int RANK>
ZVector<ValueT,RANK>::
    ZVector(std::vector<char> &_zbuf)
        : zbuf(_zbuf),
        is(std::ios_base::out | std::ios_base::binary),
        zis(is)
    {
        is.swap_vector(zbuf);

        // Read the algo from the stream
        int ialgo;
        zis.read((char *)&ialgo, sizeof(ialgo));
        boost::endian::big_to_native_inplace(ialgo);
        algo = (ZVAlgo)ialgo;

        for (int i=0; i<RANK; ++i) cur_raws[i] = 0;
    }


template<class ValueT, int RANK>
bool ZVector<ValueT,RANK>::
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
ZVector<ValueT,RANK>::
    ~ZVector()
    {
        // Get our original vector back!
        is.swap_vector(zbuf);
    }

}

}    // Namespace spsparse
#endif    // guard
