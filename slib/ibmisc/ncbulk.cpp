#include <ibmisc/ncbulk.hpp>

namespace ibmisc {

/** @param _vars {varname=interal name, fname=NetCDF filename, vname = NetcDF variable name, ...} */
NcBulkReader::NcBulkReader(
    FileLocator const *_files,
    std::vector<std::string> const &_vars)
: files(_files)
{
    size_t n = _vars.size();
    if (3*(n/3) != n) (*ibmisc_error)(-1,
        "NcBulkReader initializers must be in triplets: <varname>, <fname>, <vname>");

    for (auto ii = _vars.begin(); ii != _vars.end(); ) {
        std::string const &var(*ii++);
        std::string const &fname(*ii++);
        std::string const &vname(*ii++);

        varmap.insert(std::make_pair(var, std::array<std::string,2>{fname, vname}));
    }
}

void NcBulkReader::operator()()
{
    // Check that every expected variable has been assigned.
    if (varmap.size() != 0) {
        for (auto ii=varmap.begin(); ii != varmap.end(); ++ii)
            fprintf(stderr, "    Unassigned: %s\n", ii->first.c_str());
        (*ibmisc_error)(-1, "Unassigned variables");
    }


    // Read variables, grouped by file
    std::sort(actions.begin(), actions.end());
    std::unique_ptr<NcIO> ncio;
    std::string last_fname = "";
    for (auto ii=actions.begin(); ii != actions.end(); ++ii) {
        if (ii->fname != last_fname) {
            std::string located_fname(files->locate(ii->fname));
//printf("Opening %s -> %s\n", ii->fname.c_str(), located_fname.c_str());
            ncio.reset(new NcIO(located_fname, 'r'));
        }
        ii->fn(*ncio);
        last_fname = ii->fname;
    }
}

}    // namespace ibmisc
