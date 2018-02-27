#ifndef IBMISC_RLARRAY_HPP
#define IBMISC_RLARRAY_HPP

#include <cmath>
#include <blitz/array.h>
#include <spsparse/vector.hpp>
#include <spsparse/runlength.hpp>
#include <spsparse/sparsearray.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/stdio.hpp>

namespace ibmisc {

// --------------------------------------------------------------------------
/** Data structure for the storage of ONE runlength-encoded vector
NOTE: Count and Value here are singluar because this is the inner-most data structure. */
template<class CountT, class ValueT, class EqualT=spsparse::DefaultRLEqual<ValueT>>
struct RLVector {
    std::vector<CountT> count;
    std::vector<ValueT> value;
    spsparse::RLAlgo algo;    // Which variant of runlength encoding this is...
    EqualT eq;        // Compare when encoding

    RLVector(EqualT const &&_eq=EqualT()) : eq(std::move(_eq)) {}

    RLVector(spsparse::RLAlgo _algo, EqualT const &&_eq=EqualT()) //DefaultRLEqual<typename ValueT>())
        : algo(_algo), eq(std::move(_eq)) {}

    void ncio(NcIO &ncio, std::string const &vname,
        std::string const &count_snc_type=get_nc_type<CountT>(),
        std::string const &value_snc_type=get_nc_type<ValueT>());

    typedef spsparse::vaccum::RLEncode<
        spsparse::vaccum::Vector<CountT>,
        spsparse::vaccum::Vector<ValueT>,
        EqualT> vaccum_type;

    /** Returns a VAccum that runlength encodes into this RLVector. */
    inline vaccum_type vaccum()
    {
        return spsparse::vaccum::rl_encode(
            spsparse::vaccum::vector(count),
            spsparse::vaccum::vector(value),
            algo, EqualT(eq));
    }

    typedef spsparse::RLDecode<
        typename std::vector<CountT>::const_iterator,
        typename std::vector<ValueT>::const_iterator>
        generator_type;

    generator_type generator() const
    {
        return generator_type(
            count.cbegin(), count.cend(),
            value.cbegin(), value.cend(),
            algo);
    }

};
// --------------------------------------------------------------------------
template<class CountT, class IndexT, class ValueT, int RANK>
class RLSparseArray;

template<class CountT, class IndexT, class ValueT, int RANK>
class RLSparseArray_Generator {
    std::vector<typename RLVector<CountT,IndexT>::generator_type> indices_ii;
    typename RLVector<CountT,ValueT>::generator_type value_ii;

    static std::vector<typename RLVector<CountT,IndexT>::generator_type> make_indices_ii(RLSparseArray<CountT,IndexT,ValueT,RANK> const &arr)
    {
        std::vector<typename RLVector<CountT,IndexT>::generator_type> ret;
        for (int i=0; i<RANK; ++i) ret.push_back(arr.indices[i].generator());
        return ret;
    }

public:
    RLSparseArray_Generator(RLSparseArray<CountT,IndexT,ValueT,RANK> const &arr);

    bool operator++();

    IndexT index(int ix) const
        { return *indices_ii[ix]; }

    std::array<IndexT,RANK> index() const;

    ValueT value() const
        { return *value_ii; }

};


/** NOTE: Singular used for template parameters, plural used for
    vectors holding plural things (or accumulators).
TODO: Add a ValueEqualT template parameter, to allow the user to change it if needed (it's not needed). */
template<class CountT, class IndexT, class ValueT, int RANK>
class RLSparseArray {

    friend class RLSparseArray_Generator<CountT,IndexT,ValueT,RANK>;

    std::vector<RLVector<CountT,IndexT>> indices;
    RLVector<CountT,ValueT> values;
    int _nnz;    // Number of elements added to this RLSparseArray
    std::array<long, RANK> _shape;


    static std::vector<RLVector<CountT,IndexT>> make_indices()
    {
        std::vector<RLVector<CountT,IndexT>> ret;
        for (int i=0; i<RANK; ++i) ret.push_back(RLVector<CountT,IndexT>(spsparse::RLAlgo::DIFFS));
        return ret;
    }

public:
    RLSparseArray(std::array<long, RANK> shape) :
        indices(make_indices()),
        values(RLVector<CountT,ValueT>(spsparse::RLAlgo::PLAIN)),
        _nnz(0), _shape(shape)
    {}


    int nnz() const { return _nnz; }
    std::array<long,RANK> const &shape() const { return _shape; }

    void ncio(NcIO &ncio, std::string const &vname,
        std::string const &count_snc_type=get_nc_type<CountT>(),
        std::string const &index_snc_type=get_nc_type<IndexT>(),
        std::string const &value_snc_type=get_nc_type<ValueT>());

    typedef spsparse::accum::SparseArrayAgg<
        typename RLVector<CountT,IndexT>::vaccum_type,
        typename RLVector<CountT,ValueT>::vaccum_type,
        RANK> accum_type;

    /** Creates an encoder */
    accum_type accum();

    typedef RLSparseArray_Generator<CountT,IndexT,ValueT,RANK> generator_type;

    generator_type generator() const
        { return generator_type(*this); }


    template<class AccumT>
    void spcopy(AccumT &&ret) const
    {
        for (auto gen(generator()); ++gen;) {
            ret.add(gen.index(), gen.value());
        }
    }

};

// ==================================================================
template<class CountT, class ValueT, class EqualT>
void RLVector<CountT,ValueT,EqualT>::
    ncio(NcIO &ncio, std::string const &vname,
        std::string const &count_snc_type,
        std::string const &value_snc_type)
    {
        auto dims(get_or_add_dims(ncio, count, {vname + ".rlsize"}));    // size ignored on read

        // No need to track and record nnz; but if we do, this is how
        auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        get_or_put_att_enum(info_v, ncio.rw, "algo", algo);

        netCDF::NcVar ncv;
        ncio_vector(ncio, count, true, vname+".count", count_snc_type, dims);

        netCDF::NcVar value_v(ncio_vector(ncio, value, true, vname+".value", value_snc_type, dims));
    }


// ==================================================================


template<class CountT, class IndexT, class ValueT, int RANK>
RLSparseArray_Generator<CountT,IndexT,ValueT,RANK>::
    RLSparseArray_Generator(RLSparseArray<CountT,IndexT,ValueT,RANK> const &arr) :
        indices_ii(make_indices_ii(arr)),
        value_ii(arr.values.generator())
//        value_ii(typename RLVector<CountT,ValueT>::generator_type())
    {}


template<class CountT, class IndexT, class ValueT, int RANK>
bool RLSparseArray_Generator<CountT,IndexT,ValueT,RANK>::
    operator++()
    {
        // Increment the iterators
        bool good = ++value_ii;
        for (int i=0; i<RANK; ++i) {
            bool this_good = ++indices_ii[i];
            if (this_good != good) (*ibmisc_error)(-1,
                "All iterators should be of same length (more debugging needed here)");
        }
        return good;
    }

template<class CountT, class IndexT, class ValueT, int RANK>
std::array<IndexT,RANK> RLSparseArray_Generator<CountT,IndexT,ValueT,RANK>::
    index() const
    {
        std::array<IndexT,RANK> ret;
        for (int i=0; i<RANK; ++i) ret[i] = index(i);
        return ret;
    }

// ===================================================================
template<class CountT, class IndexT, class ValueT, int RANK>
void RLSparseArray<CountT,IndexT,ValueT,RANK>::
    ncio(NcIO &ncio, std::string const &vname,
        std::string const &count_snc_type,
        std::string const &index_snc_type,
        std::string const &value_snc_type)
    {
        auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        get_or_put_att(info_v, ncio.rw, "nnz", "int64", &_nnz, 1);
        get_or_put_att(info_v, ncio.rw, "shape", "int64", &_shape[0], RANK);

        for (int i=0; i<RANK; ++i) {
            indices[i].ncio(ncio,
                ibmisc::strprintf("%s.indices_%d", vname.c_str(), i),
                count_snc_type, index_snc_type);
        }
        values.ncio(ncio, vname+".values",
            count_snc_type, value_snc_type);
    }


template<class CountT, class IndexT, class ValueT, int RANK>
typename RLSparseArray<CountT,IndexT,ValueT,RANK>::accum_type
RLSparseArray<CountT,IndexT,ValueT,RANK>::
    accum()
    {
        std::vector<typename RLVector<CountT,IndexT>::vaccum_type> vaccum_indices;
        for (int i=0; i<RANK; ++i) {
            vaccum_indices.push_back(indices[i].vaccum());
        }

        return accum_type(std::move(vaccum_indices),
            values.vaccum(),
            _nnz);
    }

// ==================================================================








// ---------------------------------------------------------------------
} // namespace ibmisc
#endif
