#include <cstdint>
#include <ibmisc/error.hpp>
#include <ibmisc/endian.hpp>
#include <boost/endian/conversion.hpp>

namespace ibmisc {

void big_to_native(char *buf, int const item_size, long nitem)
{
    if (item_size == 1) return;

    for (long i=0; i<nitem; ++i) {
        switch(item_size) {
            case 2:
                boost::endian::big_to_native_inplace<uint16_t>(
                    *(uint16_t *)(buf + item_size*i));
            break;
            case 4:
                boost::endian::big_to_native_inplace<uint32_t>(
                    *(uint32_t *)(buf + item_size*i));
            break;
            case 8:
                boost::endian::big_to_native_inplace<uint64_t>(
                    *(uint64_t *)(buf + item_size*i));
            break;
            default:
                (*ibmisc_error)(-1, "Illegal item size: %d", item_size);
            return;
        }
    }
}



void little_to_native(char *buf, int const item_size, long nitem)
{
    if (item_size == 1) return;

    for (long i=0; i<nitem; ++i) {
        switch(item_size) {
            case 2:
                boost::endian::little_to_native_inplace<uint16_t>(
                    *(uint16_t *)(buf + item_size*i));
            break;
            case 4:
                boost::endian::little_to_native_inplace<uint32_t>(
                    *(uint32_t *)(buf + item_size*i));
            break;
            case 8:
                boost::endian::little_to_native_inplace<uint64_t>(
                    *(uint64_t *)(buf + item_size*i));
            break;
            default:
                (*ibmisc_error)(-1, "Illegal item size: %d", item_size);
            return;
        }
    }
}

void endian_to_native(char *buf, int const item_size, long nitem, Endian endian)
{
    if (endian == Endian::BIG) {
        ibmisc::big_to_native(buf, item_size, nitem);
    } else {
        ibmisc::little_to_native(buf, item_size, nitem);
    }
}



}     // namespace

// =================================================================
namespace boost {

template<>
ibmisc::Endian lexical_cast<ibmisc::Endian, std::string>(std::string const &token)
{
    if (token == "little")
        return ibmisc::Endian::LITTLE;
    if (token == "big")
        return ibmisc::Endian::BIG;

    throw boost::bad_lexical_cast();    
}

template<>
std::string lexical_cast<std::string, ibmisc::Endian>(ibmisc::Endian const &endian)
{
    switch(endian) {
        case ibmisc::Endian::LITTLE:
            return "little";
        case ibmisc::Endian::BIG:
            return "big";
    }

    throw boost::bad_lexical_cast();    
}

}    // namespace boost
// =================================================================
