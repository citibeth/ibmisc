#ifndef IBMISC_BUNDLE_HPP
#define IBMISC_BUNDLE_HPP

#include <ibmisc/blitz.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/IndexSet.hpp>

namespace ibmisc {


// ===============================================================
/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
template<class TypeT, int RANK>
class ArrayBundle {
    // Stores memory for arrays allocated as a multi-array
    TmpAlloc tmp;
public:

    struct Data {
        ArrayMeta<RANK> meta;
        blitz::Array<TypeT,RANK> arr;

#if 0
        Data(ArrayMeta<RANK> &&meta,
            blitz::Array<TypeT,RANK> &arr)
        : meta(std::move(_meta)), arr(arr) {}
#endif

        Data(
            std::string const &_name,
            blitz::Array<TypeT,RANK> const &_arr,
            std::array<int, RANK> const &_shape,
            std::array<std::string,RANK> _sdims,
            std::vector<std::tuple<std::string, std::string>> &&_attr)
        : meta(ArrayMeta<RANK>(_name,_shape, _sdims, std::move(_attr))),
            arr(_arr) {}


        void allocate(bool check,
            blitz::GeneralArrayStorage<RANK> const &storage)
        {
            if (check && arr.data()) (*ibmisc_error)(-1,
                "ArrayBundle variable %s already allocated", meta.name.c_str());
            arr.reference(blitz::Array<TypeT,RANK>(ibmisc::to_tiny<int,int,RANK>(meta.shape), storage));
        }

        std::vector<NamedDim> named_dims()
        {
            std::vector<NamedDim> ret;
            for (int i=0; i<RANK; ++i) {
                ret.push_back(NamedDim(meta.name, meta.shape));
            }
            return ret;
        }
    };

    // -------------------------------------------------------------------

    IndexSet<std::string> index;
    std::vector<Data> data;

    blitz::Array<TypeT, RANK> const &array(std::string const &name) const
        { return data[index.at(name)].arr; }

    blitz::Array<TypeT, RANK> &array(std::string const &name)
        { return data[index.at(name)].arr; }


    Data const &at(std::string const &name) const
        { return data[index.at(name)]; }

    Data &at(std::string const &name)
        { return data[index.at(name)]; }

private:

    static std::vector<std::tuple<std::string, std::string>> make_attrs(
        std::initializer_list<std::string> const &vattr);

public:
    static Data def(
        std::string const &name,
        std::initializer_list<std::string> const &vattr);

    static Data def(
        std::string const &name,
        std::array<int, RANK> const &shape,
        std::array<std::string,RANK> sdims,
        std::initializer_list<std::string> const &vattr);

    ArrayBundle() {}

    ArrayBundle(std::vector<Data> _data);


    void add(
        std::string const &name,
        std::initializer_list<std::string> const &vattr);

    void add(
        std::string const &name,
        std::array<int, RANK> const &shape,
        std::array<std::string,RANK> sdims,
        std::initializer_list<std::string> const &vattr);

    void add(
        std::string const &name,
        blitz::Array<TypeT, RANK> &arr,
        std::array<std::string,RANK> sdims,
        std::initializer_list<std::string> const &vattr);


    // ------------------------------------------------------------------
    // Allocate All Variables in a Bundle

    void set_shape(
        std::array<int, RANK> const &shape,
        std::array<std::string,RANK> sdims,
        bool check = true);

    void allocate(bool check = true,
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

    void allocate(
        std::array<int, RANK> const &_shape,
        std::array<std::string,RANK> sdims,
        bool check = true,
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

    // ------------------------------------------------------------------
    // Allocate Some Variables in a Bundle

    void set_shape(
        std::vector<std::string> const &vnames,
        std::array<int, RANK> const &shape,
        std::array<std::string,RANK> sdims,
        bool check = true);

    void allocate(
        std::vector<std::string> const &vnames,
        bool check = true,
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

    void allocate(
        std::vector<std::string> const &vnames,
        std::array<int, RANK> const &_shape,
        std::array<std::string,RANK> sdims,
        bool check = true,
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

    // -------------------------------------------------------------------
private:
#   define NCIO_BUNDLE_PARAMS \
        NcIO &ncio, \
        std::vector<std::string> const &vars, \
        std::string const &prefix, \
        std::string const &snc_type
#   define NCIO_BUNDLE_ARGS ncio, vars, prefix, snc_type

    typedef std::function<netCDF::NcVar (NCIO_BLITZ_PARAMS, std::vector<std::string> const &)> ncio_blitz_fn;

    void ncio(
        NCIO_BUNDLE_PARAMS,
        ncio_blitz_fn const &_ncio_blitz_fn);

public:
    void ncio(
        NCIO_BUNDLE_PARAMS,
        std::vector<netCDF::NcDim> const &ncdims={},
        DimOrderMatch match=DimOrderMatch::MEMORY,
        bool ncdims_in_nc_order=true);

    void ncio_alloc(
        NCIO_BUNDLE_PARAMS,
        std::vector<netCDF::NcDim> const &ncdims = {},
        DimOrderMatch match=DimOrderMatch::MEMORY,
        bool ncdims_in_nc_order = true,
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>());

    void ncio_alloc(
        NCIO_BUNDLE_PARAMS,
        DimOrderMatch match=DimOrderMatch::MEMORY,
        blitz::GeneralArrayStorage<RANK> const &storage = blitz::GeneralArrayStorage<RANK>())
    {
        ncio_alloc(NCIO_BUNDLE_ARGS, {}, match, true, storage);
    }

    void ncio_partial(
        NCIO_BUNDLE_PARAMS,
        std::vector<netCDF::NcDim> const &ncdims,
        std::vector<int> const &nc_start,    // Where to start each dimension in NetCDF
        std::array<int,RANK> const &b2n);    // Where to slot each Blitz++ dimension



};
// --------------------------------------------------------------------
template<class TypeT, int RANK>
std::vector<std::tuple<std::string, std::string>> ArrayBundle<TypeT,RANK>::make_attrs(
    std::initializer_list<std::string> const &vattr)
{
    std::vector<std::tuple<std::string, std::string>> ret;
    for (auto ii = vattr.begin(); ii != vattr.end(); ) {
        std::string const &key(*ii++);
        if (ii == vattr.end()) (*ibmisc_error)(-1,
            "Odd number of strings in (key,value) attr list");
        std::string const &value(*ii++);

        ret.push_back(std::make_tuple(key, value));
    }
    return ret;
}

template<class TypeT, int RANK>
typename ArrayBundle<TypeT,RANK>::Data ArrayBundle<TypeT,RANK>::def(
    std::string const &name,
    std::initializer_list<std::string> const &vattr)
{
    std::array<int,RANK> shape;
    std::array<std::string,RANK> sdims;
    for (int i=0; i<RANK; ++i) shape[i] = -1;
    return Data(name, blitz::Array<TypeT,RANK>(), shape,
        std::move(sdims), make_attrs(vattr));
}

template<class TypeT, int RANK>
typename ArrayBundle<TypeT,RANK>::Data ArrayBundle<TypeT,RANK>::def(
    std::string const &name,
    std::array<int, RANK> const &shape,
    std::array<std::string,RANK> sdims,
    std::initializer_list<std::string> const &vattr)
{
    return Data(name, blitz::Array<TypeT,RANK>(), shape,
        std::move(sdims), make_attrs(vattr));
}

template<class TypeT, int RANK>
ArrayBundle<TypeT,RANK>::ArrayBundle(std::vector<Data> _data) : data(std::move(_data))
{
    for (Data &meta : data) {
        index.insert(meta.meta.name);
    }
}




template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    std::initializer_list<std::string> const &vattr)
{
    data.push_back(def(name, vattr));
    index.insert(data.back().meta.name);
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    std::array<int, RANK> const &shape,
    std::array<std::string,RANK> sdims,
    std::initializer_list<std::string> const &vattr)
{
    data.push_back(def(name, shape, std::move(sdims), vattr));
    index.insert(data.back().meta.name);
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::Array<TypeT,RANK> &arr,
    std::array<std::string,RANK> sdims,
    std::initializer_list<std::string> const &vattr)
{
    data.push_back(Data(name, arr, arr.shape(),
        std::move(sdims), make_attrs(vattr)));
}


// ------------------------------------------------------------------
// Allocate All Variables in a Bundle

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::set_shape(
    std::array<int, RANK> const &shape,
    std::array<std::string,RANK> sdims,
    bool check)
{
    for (auto &meta : data) {
        if (meta.shape[0] < 0)
            meta.meta.set_shape(shape, sdims, check);
    }
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::allocate(bool check,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    for (auto &meta : data) {
        if (!meta.arr.data())
            meta.allocate(check, storage);
    }
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::allocate(
    std::array<int, RANK> const &_shape,
    std::array<std::string,RANK> sdims,
    bool check,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    for (auto &meta : data) {
        if (meta.meta.shape[0] < 0 || !meta.arr.data()) {
            meta.meta.set_shape(_shape, sdims, check);
            meta.allocate(check, storage);
        }
    }
}

// ------------------------------------------------------------------
// Allocate Some Variables in a Bundle

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::set_shape(
    std::vector<std::string> const &vnames,
    std::array<int, RANK> const &shape,
    std::array<std::string,RANK> sdims,
    bool check)
{
    for (auto &vname : vnames) {
        auto &meta(data[index.at(vname)]);
        meta.meta.set_shape(shape, sdims, check);
    }
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::allocate(
    std::vector<std::string> const &vnames,
    bool check,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    for (auto &vname : vnames) {
        auto &meta(data[index.at(vname)]);
        meta.allocate(check, storage);
    }
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::allocate(
    std::vector<std::string> const &vnames,
    std::array<int, RANK> const &_shape,
    std::array<std::string,RANK> sdims,
    bool check,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    for (auto &vname : vnames) {
        auto &meta(data[index.at(vname)]);
        meta.meta.set_shape(_shape, sdims, check);
        meta.allocate(check, storage);
    }
}
// --------------------------------------------------------------------
// --------------------------------------------------------------------
template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio(
    NCIO_BUNDLE_PARAMS,
    ncio_blitz_fn const &_ncio_blitz_fn)
{
    // Determine which variables in the bundle we will operate on
    std::vector<std::string> all_vars;
    std::vector<std::string> const *myvars;
    if (vars.size() == 0) {
        for (size_t i=0; i<index.size(); ++i) {
            all_vars.push_back(index[i]);
        }
        myvars = &all_vars;
    } else {
        myvars = &vars;
    }


    // Do ncio_blitz() for each varaible
    for (auto &var : *myvars) {
        int i=index.at(var);

        auto &meta(data[i]);
        std::string vname(prefix+meta.meta.name);

        // Delegate to lower level to read/write this array.
        auto &arr(meta.arr);    // Used by NCIO_BLITZ_ARGS macro
        _ncio_blitz_fn(NCIO_BLITZ_ARGS, to_vector(meta.meta.sdims));

        // Read/write attributes
        netCDF::NcVar ncvar = ncio.nc->getVar(vname);
        if (ncio.rw == 'w') {
            for (auto &kv : meta.meta.attr) {
                std::string const &name(std::get<0>(kv));
                std::string &value(std::get<1>(kv));    // Read back into Data

                get_or_put_att(ncvar, ncio.rw, name, value);
            }
        } else {
            meta.meta.attr.clear();
            auto atts(ncvar.getAtts());
            for (auto ii(atts.begin()); ii != atts.end(); ++ii) {
                std::string const &aname(ii->first);
                std::string aval;
                ii->second.getValues(aval);

                meta.meta.attr.push_back(std::make_tuple(aname, aval));
            }
        }
    }
}
// --------------------------------------------------------------------
template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio(
    NCIO_BUNDLE_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    DimOrderMatch match,
    bool ncdims_in_nc_order)
{
    using namespace std::placeholders;
    auto fn(std::bind(&ncio_blitz<TypeT,RANK>, _1, _2, _3, _4,
            ncdims, match, ncdims_in_nc_order, _5));
    this->ncio(NCIO_BUNDLE_ARGS, fn);
}

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio_alloc(
    NCIO_BUNDLE_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    DimOrderMatch match,
    bool ncdims_in_nc_order,
    blitz::GeneralArrayStorage<RANK> const &storage)
{
    using namespace std::placeholders;
    this->ncio(NCIO_BUNDLE_ARGS,
        std::bind(&ncio_blitz_alloc<TypeT,RANK>, _1, _2, _3, _4,
            ncdims, match, ncdims_in_nc_order, storage, _5));
}


template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio_partial(
    NCIO_BUNDLE_PARAMS,
    std::vector<netCDF::NcDim> const &ncdims,
    std::vector<int> const &nc_start,    // Where to start each dimension in NetCDF
    std::array<int,RANK> const &b2n)    // Where to slot each Blitz++ dimension
{
    using namespace std::placeholders;
    this->ncio(NCIO_BUNDLE_ARGS,
        std::bind(&ncio_blitz_partial<TypeT,RANK>, _1, _2, _3, _4,
            ncdims, nc_start, b2n, _5));
}


// -------------------------------------------------------------
/** Reshapes a bundle of Blitz++ arrays to a bundle of 1-D Blitz++ array */
template<class TypeT, int RANK>
ArrayBundle<TypeT,1> reshape1(
    ArrayBundle<TypeT, RANK> &bundle,
    int lbound = 0,
    std::array<std::string,1> const &sdims = {""});

template<class TypeT, int RANK>
ArrayBundle<TypeT,1> reshape1(
    ArrayBundle<TypeT, RANK> &bundle,
    int lbound = 0,
    std::array<std::string,1> const &sdims = {""})
{
    ArrayBundle<TypeT,1> bundle1;
    for (size_t i=0; i<bundle.index.size(); ++i) {
        auto &meta(bundle.data[i]);
        blitz::Array<TypeT,1> meta1_arr(reshape1(meta.arr, lbound));
        typename ArrayBundle<TypeT,1>::Data meta1(
            meta.meta.name,
            meta1_arr,
            blitz::shape(meta1_arr.extent(0)),
            sdims,
            meta.meta.attr);
        bundle1.index.insert(meta.meta.name);
        bundle1.data.push_back(meta1);
    }
    return bundle1;
}

}

#endif // IBMISC_BUNDLE_HPP
