#ifndef IBMISC_ZSPARSEARRAY_HPP
#define IBMISC_ZSPARSEARRAY_HPP

#include <cmath>
#include <blitz/array.h>
#include <spsparse/zvector.hpp>
#include <spsparse/sparsearray.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/stdio.hpp>

namespace ibmisc {

// ======================================================================

/** Output is flushed once the ZArray_Accum object is destroyed. */
template<class IndexT, class ValueT, int RANK>
class ZArray_Accum {
public:
    typedef ValueT val_type;
    typedef IndexT index_type;
    static int const rank = RANK;

private:
    spsparse::vaccum::ZVector<IndexT,RANK> indices;
    spsparse::vaccum::ZVector<ValueT,1> values;
    std::array<long,RANK> &shape;
    long &nnz;

public:
    ZArray_Accum(
        std::vector<char> &_indices,    // Holds std::array<IndexT,RANK>
        std::vector<char> &_values,     // Holds std::array<ValueT,1>
        std::array<long,RANK> &_shape,
        long &_nnz);

    void add(std::array<IndexT,RANK> const &index, ValueT const &value)
    {
        ++nnz;
        indices.add(index);
        values.add({value});
    }

    void set_shape(std::array<long,RANK> const &_shape)
        { shape = _shape; }
};


template<class IndexT, class ValueT, int RANK>
ZArray_Accum<IndexT,ValueT,RANK>::
    ZArray_Accum(
        std::vector<char> &_indices,    // Holds std::array<IndexT,RANK>
        std::vector<char> &_values,     // Holds std::array<ValueT,1>
        std::array<long,RANK> &_shape,
        long &_nnz)
    : indices(_indices, spsparse::ZVAlgo::DIFFS),
        values(_values, spsparse::ZVAlgo::PLAIN),
        shape(_shape), nnz(_nnz)
    {}


// ========================================================================

template<class IndexT, class ValueT, int RANK>
class ZArray_Generator {
    spsparse::vgen::ZVector<IndexT,RANK> indices;
    spsparse::vgen::ZVector<ValueT,1> values;

public:
    ZArray_Generator(
        std::vector<char> &_indices,    // Holds std::array<IndexT,RANK>
        std::vector<char> &_values)     // Holds std::array<ValueT,1>
    : indices(_indices), values(_values) {}

    bool operator++();

    IndexT index(int ix) const
        { return (*indices)[ix]; }

    std::array<IndexT,RANK> const &index() const
        { return *indices; }

    ValueT value() const
        { return (*values)[0]; }


    ZArray_Generator<IndexT,ValueT,RANK> const *operator->() const { return this; }
};

template<class IndexT, class ValueT, int RANK>
bool ZArray_Generator<IndexT,ValueT,RANK>::
    operator++()
    {
        // Increment the iterators
        bool good_indices = ++indices;
        bool good_values = ++values;
        if (good_indices != good_values) (*ibmisc_error)(-1,
            "All iterators should be of same length (more debugging needed here)");
		// if (good_indices) printf("   ++  (%d, %d) %g\n", (*indices)[0], (*indices)[1], (*values)[0]);
        return good_indices;
    }


// =========================================================================

/** A compressed sparse array format.
Most directly, a compressed verison of TupleList.

NOTE: Singular used for template parameters, plural used for
    vectors holding plural things (or accumulators).
TODO: Add a ValueEqualT template parameter, to allow the user to change it if needed (it's not needed). */
template<class IndexT, class ValueT, int RANK>
class ZArray {

    friend class ZArray_Generator<IndexT,ValueT,RANK>;

    std::vector<char> indices;    // Holds std::array<IndexT,RANK>
    std::vector<char> values;     // Holds std::array<ValueT,1>
    long _nnz;    // Number of elements added to this ZArray
    std::array<long, RANK> _shape;

public:
    ZArray() : _nnz(0)
    {
        for (int i=0; i<RANK; ++i) _shape[i] = 0;
    }

    ZArray(std::array<long, RANK> const &shape) :
        _nnz(0), _shape(shape)
    {}

    long nnz() const { return _nnz; }
    std::array<long,RANK> const &shape() const { return _shape; }

    void ncio(NcIO &ncio, std::string const &vname);

    typedef ZArray_Accum<IndexT,ValueT,RANK> accum_type;
    /** Creates an encoder */
    accum_type accum()
        { return accum_type(indices, values, _shape, _nnz); }


    typedef ZArray_Generator<IndexT,ValueT,RANK> generator_type;
    generator_type generator() const
    {
        return generator_type(
            *const_cast<std::vector<char> *>(&indices),
            *const_cast<std::vector<char> *>(&values));
    }


    template<class AccumT>
    void spcopy(AccumT &&ret) const
    {
        for (auto gen(generator()); ++gen;) {
            ret.add(gen.index(), gen.value());
        }
    }

};


template<class IndexT, class ValueT, int RANK>
void ZArray<IndexT,ValueT,RANK>::
    ncio(NcIO &ncio, std::string const &vname)
    {
        auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        std::string zarray("ZARRAY");
        get_or_put_att(info_v, ncio.rw, "format", zarray);
        int rank = RANK;
        get_or_put_att(info_v, ncio.rw, "rank", "int", &rank, 1);
        get_or_put_att(info_v, ncio.rw, "nnz", "int64", &_nnz, 1);
        get_or_put_att(info_v, ncio.rw, "shape", "int64", &_shape[0], RANK);

        netCDF::NcVar ncvar;
        ncvar = ncio_vector<char,uint8_t>(
            ncio, indices, true, vname+".indices", "ubyte",
            get_or_add_dims(ncio, indices, {vname + ".indices.zsize"}));
        if (ncio.rw == 'w')
            ncvar.setCompression(false, false, 0);    // We're already compressing, NetCDF should not also compress

        ncvar = ncio_vector<char,uint8_t>(
            ncio, values, true, vname+".values", "ubyte",
            get_or_add_dims(ncio, values, {vname + ".values.zsize"}));
        if (ncio.rw == 'w')
            ncvar.setCompression(false, false, 0);    // We're already compressing, NetCDF should not also compress
    }

// ---------------------------------------------------------------------
} // namespace ibmisc
#endif
