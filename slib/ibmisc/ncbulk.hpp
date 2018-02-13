#ifndef IBMISC_NCBULK_HPP
#define IBMISC_NCBULK_HPP

#include <map>
#include <functional>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/filesystem.hpp>

namespace ibmisc {


/** Reads a bunch of blitz::Arrays from multiple NetCDF files.  It is
    flexible about which arrays come from which files.  Specs are
    provided as follows:


    1. The end user provides a list of string tuples (arrayname,
       filename, ncvarname).  This specifies that the array named
       `arrayname` should be loaded from NetCDF variable `ncvarname`
       file `filename`.  This list is provided to the constructor.

    2. The programmer calls operator()(arrayname, blitz::Array)
       repeatedly, once per blitz::Array that needs to be loaded.

    3. The programmer calls operator()(), or the destructur.
       NcBulkReader then matches up the specs provided by the operator
       and end user, opens up each files once, and loads the specified
       arrays from them.
*/
class NcBulkReader
{
    FileLocator const * const files;
    std::map<std::string, std::array<std::string,2>> varmap;

    // Actions, keyed by filename
    typedef std::function<netCDF::NcVar (NcIO &ncio)> ActionFn;
    struct Action {
        std::string fname;
        ActionFn fn;

        Action(std::string const &_fname, ActionFn const &_fn)
            : fname(_fname), fn(_fn) {}

        bool operator<(Action const &other) const
            { return fname < other.fname; }
    };
    std::vector<Action> actions;

public:

    /** @param _vars {varname=interal name, fname=NetCDF filename, vname = NetcDF variable name, ...} */
    NcBulkReader(
        FileLocator const *_files,
        std::vector<std::string> const &_vars);

private:

    template<class TypeT, int RANK>
    void _add_var(
        blitz::Array<TypeT, RANK> &var,
        std::string const &fname,
        std::string const &vname);


public:

    /** Add to the set of Blitz variables to be read.
    @param varname Internal name for this variable, assigned in
        constructor */
    template<class TypeT, int RANK>
    NcBulkReader &operator()(
        std::string const &varname,
        blitz::Array<TypeT, RANK> &var);

    /** Do the read! */
    void operator()();
    ~NcBulkReader() { this->operator()(); }
};

// ----------------------------------------------------
template<class TypeT, int RANK>
void NcBulkReader::_add_var(
    blitz::Array<TypeT, RANK> &var,
    std::string const &fname,
    std::string const &vname)
{
    ActionFn action(std::bind(
        &ncio_blitz<TypeT,RANK>, std::placeholders::_1, var, vname,
        "", std::vector<netCDF::NcDim>{},
        DimOrderMatch::MEMORY, true, std::vector<std::string>{}
    ));

    actions.push_back(Action(fname, action));

}

template<class TypeT, int RANK>
NcBulkReader &NcBulkReader::operator()(
    std::string const &varname,
    blitz::Array<TypeT, RANK> &var)
{
    auto iix(varmap.find(varname));
    if (iix == varmap.end()) (*ibmisc_error)(-1,
        "User failed to include fname and vname for variable %s", varname.c_str());

    std::string const &fname(iix->second[0]);
    std::string const &vname(iix->second[1]);

    _add_var(var, fname, vname);
    varmap.erase(iix);    // We've added once, can't add again

    return *this;
}


}    // namespace ibmisc
#endif // IBMISC_NCBULK_HPP
