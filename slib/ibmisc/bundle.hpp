#ifndef IBMISC_BUNDLE_HPP
#define IBMISC_BUNDLE_HPP

#include <ibmisc/blitz.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/IndexSet.hpp>

namespace ibmisc {


/** Array meta-data stored with GISS-type files */
template<class TypeT, int RANK>
struct ArrayMeta {
    std::string const name;
    blitz::Array<TypeT, RANK> arr;
    std::array<std::string,RANK> const sdims;
    std::vector<std::tuple<std::string, std::string>> attr;    // (attr-name, value)

    ArrayMeta(
        std::string const &_name,
        blitz::Array<TypeT, RANK> const &_arr,
        std::array<std::string,RANK> const &_sdims,
    std::vector<std::tuple<std::string, std::string>> const &_attr)
    : name(_name), arr(_arr), sdims(_sdims), attr(_attr) {}
};
// -----------------------------------------------------------------
struct BundleSpec {
    std::string name;
    std::vector<std::tuple<std::string, std::string>> attr;
    BundleSpec(std::string const &_name,
        std::vector<std::tuple<std::string, std::string>> const &_attr)
    : name(_name), attr(_attr) {}
};

/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
template<class TypeT, int RANK>
class ArrayBundle {
    // Stores memory for arrays allocated as a multi-array
    std::vector<blitz::Array<TypeT, RANK+1>> multi_arr;
public:
    IndexSet<std::string> index;
    std::vector<ArrayMeta<TypeT, RANK>> data;

public:
    blitz::Array<TypeT, RANK> const &at(std::string const &name) const
        { return data[index.at(name)].arr; }

    blitz::Array<TypeT, RANK> &at(std::string const &name)
        { return data[index.at(name)].arr; }


    std::vector<std::string> const &keys() const
        { return index.keys(); }

    /** Allocate an array and add */
    blitz::Array<TypeT,RANK> add(
        std::string const &name,
        blitz::TinyVector<int,RANK> const &shape,
        std::array<std::string,RANK> const &sdims,
        std::vector<std::tuple<std::string, std::string>> const &attr,    // name,value
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

    /** Add an existing array */
    blitz::Array<TypeT,RANK> add(
        std::string const &name,
        blitz::Array<TypeT, RANK> arr,
        std::array<std::string,RANK> const &sdims,
        std::vector<std::tuple<std::string, std::string>> const &attr);

    /** Allocate and add many as part of a RANK+1 array */
    blitz::Array<TypeT,RANK+1> add(
        blitz::TinyVector<int,RANK> const &shape,
        std::array<std::string,RANK> const &sdims,
        std::vector<BundleSpec> const &inits);

#if 0
    int add(ArrayMeta<TypeT, RANK> const &datum)
    {
        size_t ix = index.insert(datum.name);
        data.push_back(datum);
    }
#endif

    void ncio(NcIO &ncio, std::string const &prefix, std::string const &snc_type, std::vector<std::string> const &vars = {"<all>"});

};
// -----------------------------------------------------------------

/** Add a self-allocated array */
template<class TypeT, int RANK>
blitz::Array<TypeT,RANK> ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::TinyVector<int, RANK> const &shape,
    std::array<std::string,RANK> const &sdims,
    std::vector<std::tuple<std::string, std::string>> const &vattr,    // name,value
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    size_t ix = index.insert(name);
    // No array provided, allocate a new one
    blitz::Array<TypeT,RANK> arr(shape, storage);
    data.push_back(ArrayMeta<TypeT,RANK>(
        name, arr, sdims, vattr));
    return arr;
}

/** Add an existing array */
template<class TypeT, int RANK>
blitz::Array<TypeT,RANK> ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::Array<TypeT, RANK> arr,
    std::array<std::string,RANK> const &sdims,
    std::vector<std::tuple<std::string, std::string>> const &vattr)    // name,value
{
    size_t ix = index.insert(name);
    // Array provided, reference it
    data.push_back(ArrayMeta<TypeT,RANK>(
        name, arr, sdims, vattr));
    return arr;
}

/** Allocate and add many as part of a RANK+1 array.
Generates C-arrays only (base 0, row major).  Use c_to_f() if you
want to initialize Fortran arrays */
template<class TypeT, int RANK>
blitz::Array<TypeT,RANK+1> ArrayBundle<TypeT,RANK>::add(
    blitz::TinyVector<int,RANK> const &shape,
    std::array<std::string,RANK> const &sdims,
    std::vector<BundleSpec> const &inits)
{
    // Initialize the multi-array
    int const nvar = inits.size();
    blitz::TinyVector<int,RANK+1> SHAPE;
    SHAPE[0] = nvar;
    for (int i=0; i<RANK; ++i) SHAPE[i+1] = shape[i];
    multi_arr.push_back(blitz::Array<TypeT,RANK+1>(SHAPE));
    auto &ARR(multi_arr.at(multi_arr.size()-1));

    blitz::TinyVector<int,RANK+1> LOWER;
    for (int i=0; i<RANK; ++i) LOWER[i] = ARR.lbound(i+1);

    for (int ix=0; ix<inits.size(); ++ix) {
        auto const &init(inits[ix]);

        LOWER[0] = ARR.lbound(0) + ix;
        blitz::Array<TypeT,RANK> arr(
            &ARR(LOWER), shape, blitz::neverDeleteData);

        add(init.name, arr, sdims, init.attr);
    }

    return ARR;
}


// --------------------------------------------------------------------

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio(NcIO &ncio, std::string const &prefix, std::string const &snc_type, std::vector<std::string> const &vars)
{
    std::vector<std::string> all_vars;
    std::vector<std::string> const *myvars;
    if (vars.size() == 1 && vars[0] == "<all>") {
        for (size_t i=0; i<index.size(); ++i) {
            all_vars.push_back(index[i]);
        }
        myvars = &all_vars;
    } else {
        myvars = &vars;
    }

    for (auto &var : *myvars) {
        int i=index.at(var);

        auto &meta(data[i]);

        // Set up the dimensions
        auto dims_f(get_or_add_dims(ncio,
            meta.arr,
            to_vector(meta.sdims)));

        // Read/Write the NetCDF variable
        // (will auto-reverse dims if it detects column major)
        auto ncvar(ncio_blitz(ncio, meta.arr, false, prefix + meta.name, snc_type, dims_f));

        // Read/write attributes
        for (auto &kv : meta.attr) {
            std::string const &name(std::get<0>(kv));
            std::string &value(std::get<1>(kv));    // Read back into ArrayMeta

            get_or_put_att(ncvar, ncio.rw, name, value);
        }
    }
}



}

#endif // IBMISC_BUNDLE_HPP
