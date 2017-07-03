// http://www.eecs.umich.edu/courses/eecs380/HANDOUTS/cppBinaryFileIO-2.html
#include <ibmisc/fortranio.hpp>
#include <ibmisc/endian.hpp>
#include <boost/endian/conversion.hpp>
#include <ibmisc/endian.hpp>

namespace ibmisc {
namespace fortran {

EndR endr;


SimpleBufSpec::SimpleBufSpec(char *_buf, int _item_size, long _nitem)
    : BufSpec(_item_size*_nitem), item_size(_item_size), nitem(_nitem), buf(_buf) {}
void SimpleBufSpec::read(UnformattedInput &infile)
{
    // Read into the buffer
    size_t len = item_size * nitem;
    infile.fin.read(buf, len);

     // Swap for endian
    endian_to_native(buf, item_size, nitem, infile.endian);
}




//template<class IStreamT>
void read::operator>>(EndR const &endr)
{
    // Read header
    uint32_t len0;
    infile->fin.read((char *)&len0, 4);
    if (infile->endian == ibmisc::Endian::BIG) {
        boost::endian::big_to_native_inplace(len0);
    } else {
        boost::endian::little_to_native_inplace(len0);
    }
    if (infile->fin.fail()) (*ibmisc_error)(-1,
        "Error reading len0\n");

    // Make sure desired body is right size
    if (specs.size() > 0) {
        size_t total = 0;
        for (auto &spec : specs) total += spec->len;
        if (total != len0) (*ibmisc_error)(-1,
            "Trying to read record of size %ld with pattern of size %ld",
            (long)len0, (long)total);
        // Read the bodies
        for (auto &spec : specs) {
            spec->read(*infile);
            if (infile->fin.fail()) (*ibmisc_error)(-1,
                "Error reading into a buffer");
        }
    } else {
        // Skip past this record
        infile->fin.seekg(len0, infile->fin.cur);
    }

    // Read trailer
    uint32_t len1;
    infile->fin.read((char *)&len1, 4);
    if (infile->endian == ibmisc::Endian::BIG) {
        boost::endian::big_to_native_inplace(len1);
    } else {
        boost::endian::little_to_native_inplace(len1);
    }
    if (infile->fin.fail()) (*ibmisc_error)(-1,
        "Error reading len1\n");
    if (len0 != len1) (*ibmisc_error)(-1,
        "Record len0=%d does not match len1=%d", len0, len1);
}


}}
