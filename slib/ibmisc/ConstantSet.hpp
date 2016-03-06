#pragma once

#include <ibmisc/IndexSet.hpp>
#include <ibmisc/udunits2.hpp>
#include <ibmisc/netcdf.hpp>

namespace ibmisc {

static double const nan = std::numeric_limits<double>::quiet_NaN();

class ConstantSet
{
public:
    struct Data {
        /** The "short" name of the variable */
        std::string const name;
        /** The units of the variable, in UDUNITS format. */
        std::string const units;            //!< UDUnits-compatible string
        /** A textual description of the variable, also called the "long name" */
        std::string const description;

        double val;

        Data(std::string const &_name,
            std::string const &_units,
            std::string const &_description)
        : name(_name),
        units(_units),
        description(_description),
        val(nan)
        {}
    };

    IndexSet<std::string> index;    // Densely ordered set of constant names
    std::vector<Data> data; // Meta-data and value data

protected:
    int add(
        std::string const &name, std::string const &units,
        std::string const &description);

    UTSystem const * ut_system; //!< Unit system to use for conversions

public:
    void init(UTSystem const *_ut_system)
    {
        ut_system = _ut_system;
    }

    /** @return Index of the new constant. */
    int set(
        std::string const &name,
        double val,
        std::string const &units,
        std::string const &description);

    /** Copy a constant's value and units to another constant
    @return Index of the new constant. */
    int copy(std::string const &dst_name,
        std::string const &src_name,
        std::string const &description);

#if 0
    // ---------------------------------------------
    typedef IndexSet<std::string>::iterator iterator;
    iterator begin() { return index.begin(); }
    iterator end() { return index.end(); }
    typedef IndexSet<std::string>::const_iterator const_iterator;
    const_iterator begin() const { return index.begin(); }
    const_iterator end() const { return index.end(); }
    const_iterator cbegin() const { return index.cbegin(); }
    const_iterator cend() const { return index.cend(); }
#endif

    // ---------------------------------------------
    // Reading constants
    size_t size() { return index.size(); }


    Data const &operator[](size_t ix) const
        { return data[ix]; }

    Data const &at(std::string const &name) const
        { return data[index.at(name)]; }

    double get_as(std::string const &name,
        UTUnit const &units) const;

    double get_as(std::string const &name,
        std::string const &sunits) const;



    void ncio(NcIO &ncio, std::string const &vname);
};

}

inline std::ostream &operator<<(std::ostream &out, ibmisc::ConstantSet::Data const &cf)
    { return out << "(" << cf.name << " = " << cf.val << " [" << cf.units << "])"; } 

inline std::ostream &operator<<(std::ostream &out, ibmisc::ConstantSet const &constants)
{
    for (auto ii=constants.data.begin(); ii != constants.data.end(); ++ii) {
        std::cout << *ii << std::endl;
    }
    return out;
}
