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
    size_t nbytes = item_size * nitem;
    infile.fin.read(buf, nbytes);

     // Swap for endian
    endian_to_native(buf, item_size, nitem, infile.endian);
}




//template<class IStreamT>
void read::operator>>(EndR const &endr)
{
    // Read header
    uint32_t nbytes0;
    infile->fin.read((char *)&nbytes0, 4);
    if (infile->endian == ibmisc::Endian::BIG) {
        boost::endian::big_to_native_inplace(nbytes0);
    } else {
        boost::endian::little_to_native_inplace(nbytes0);
    }
    if (infile->eof()) return;    // We're done, user must check for EOF
    if (infile->fin.fail()) (*ibmisc_error)(-1,
        "Error reading nbytes0\n");

    // Set wildcard spec sizes
    long total = 0;
    for (auto &spec : specs) {
        if (!spec->wildcard) total += spec->nbytes;
    }
    for (auto &spec : specs) {
        if (spec->wildcard) spec->set_nbytes(nbytes0 - total);
    }

    // Make sure desired body is right size
    if (specs.size() > 0) {
        size_t total = 0;
        for (auto &spec : specs) {
            total += spec->nbytes;
        }
        if (total != nbytes0) (*ibmisc_error)(-1,
            "Trying to read record of size %ld with pattern of size %ld",
            (long)nbytes0, (long)total);
        // Read the bodies
        for (auto &spec : specs) {
            spec->read(*infile);
            if (infile->fin.fail()) (*ibmisc_error)(-1,
                "Error reading into a buffer");
        }
    } else {
        // Skip past this record
        infile->fin.seekg(nbytes0, infile->fin.cur);
    }

    // Read trailer
    uint32_t nbytes1;
    infile->fin.read((char *)&nbytes1, 4);
    if (infile->endian == ibmisc::Endian::BIG) {
        boost::endian::big_to_native_inplace(nbytes1);
    } else {
        boost::endian::little_to_native_inplace(nbytes1);
    }
    if (infile->fin.fail()) (*ibmisc_error)(-1,
        "Error reading nbytes1\n");
    if (nbytes0 != nbytes1) (*ibmisc_error)(-1,
        "Record nbytes0=%d does not match nbytes1=%d", nbytes0, nbytes1);
}


}}
