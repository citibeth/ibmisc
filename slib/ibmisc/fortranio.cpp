// http://www.eecs.umich.edu/courses/eecs380/HANDOUTS/cppBinaryFileIO-2.html
#include <ibmisc/fortranio.hpp>

namespace ibmisc {
namespace fortran {

EndR endr;

//template<class IStreamT>
void read::operator>>(EndR const &endr)
{
    // Read header
    uint32_t len0;
    infile->read((char *)&len0, 4);
    if (infile->fail()) (*ibmisc_error)(-1,
        "Error reading len0\n");

    // Make sure desired body is right size
    if (specs.size() > 0) {
        size_t total = 0;
        for (BufSpec const &spec : specs) total += spec.len;
        if (total != len0) (*ibmisc_error)(-1,
            "Trying to read record of size %ld with pattern of size %ld",
            (long)len0, (long)total);

        // Read the bodies
        for (BufSpec const &spec : specs) {
            infile->read(spec.buf, spec.len);
            if (infile->fail()) (*ibmisc_error)(-1,
                "Error reading into a buffer");
        }
    } else {
        // Skip past this record
        infile->seekg(len0, infile->cur);
    }

    // Read trailer
    uint32_t len1;
    infile->read((char *)&len1, 4);
    if (infile->fail()) (*ibmisc_error)(-1,
        "Error reading len1\n");
    if (len0 != len1) (*ibmisc_error)(-1,
        "Record len0=%d does not match len1=%d", len0, len1);
}


}}
