#ifndef IBMISC_RUNLENGTH_HPP
#define IBMISC_RUNLENGTH_HPP

#include <cmath>
#include <blitz/array.h>

namespace ibmisc {

// --------------------------------------------------------------------------
/** Data structure for the storage of ONE runlength-encoded vector
NOTE: Count and Value here are singluar because this is the inner-most data structure. */
template<class CountT, class ValueT, int RANK>
struct RLVector {
    std::vector<CountT> count;
    std::vector<ValueT> value;
    RLAlgo algo;    // Which variant of runlength encoding this is...

    void ncio(NcIO &ncio, std::string const &vname, std::string const &count_snc_type, std::string const &value_snc_type);

    typedef vaccum::RLEncode<vaccum::Vector<CountsT>, vaccum::Vector<ValueT>, EqualT> vaccum_type;

    /** Returns a VAccum that runlength encodes into this RLVector. */
    template <class EqualT=std::equal_to<typename ValueT>>
    inline vaccum_type vaccum(
        RLAlgo _algo = RLAlgo::PLAIN,
        EqualT const &&eq=DefaultRLEqual<typename ValueT>())
    {
        algo = _algo;
        return ibmisc::rl_encoder(
            vaccum::vector(count),
            vaccum::vector(value),
            diff_encode, std::move(eq));
    }

    typedef RLDecode<
        typename std::vector<CountT>::iterator,
        typename std::vector<ValueT>::iterator>
        generator_type;

    generator_type generator()
    {
        return rl_decode(
            count.begin(), count.end(),
            value.begin(), value.end(),
            algo);
    }

};
// --------------------------------------------------------------------------
template<class CountT, class IndexT, class ValueT, int RANK>
class RLSparseArray;

template<class CountT, class IndexT, class ValueT, int RANK>
class RLSparseArray_Generator {
    std::vector<typename indices_rlvector_type::generator_type> indices_ii;
    typename values_rlvector_type::generator_type value_ii;

    RLSparseArray_Generator(RLSparseArray const &arr);

    bool operator++()

    IndexT index(int ix) const
        { return *indices_ii[ix]; }

    std::array<IndexT,RANK> index() const;

    ValueT value() const
        { return *values_ii; }

};


/** NOTE: Singular used for template parameters, plural used for
    vectors holding plural things (or accumulators). */
template<class CountT, class IndexT, class ValueT, int RANK>
class RLSparseArray {

    typedef RLVector<CountT,IndexT> indices_rlvector_type;
    typedef RLVector<CountT,ValueT> values_rlvector_type;

    std::array<indices_rlvector_type,RANK> indices;
    int _nnz = 0;    // Number of elements added to this RLSparseArray

public:
    int nnz() { return _nnz; }

    void ncio(NcIO &ncio, std::string const &vname,
        std::string const &count_snc_type,
        std::string const &index_snc_type,
        std::string const &value_snc_type);

    typedef accum::SparseArrayAgg<
        typename indices_rlvector_type::vaccum_type,
        typename values_rlvector_type::vaccum_type,
        RANK> accum_type;

    /** Creates an encoder */
    template<class ValueEqualT>
    accum_type accum(ValueEqualT &&value_eq);

    RLSparseArray_Generator generator() const
        { return Generator(*this); }


    template<class AccumT>
    void spcopy(AccumT &&ret) const
    {
        for (auto gen(Generator); ++gen;) {
            ret.add(gen.index(), gen.value());
        }
    }

};

// ==================================================================
template<class CountT, class ValueT, int RANK>
void RLVector<CountT,ValueT,RANK>::
    ncio(NcIO &ncio, std::string const &vname, std::string const &count_snc_type, std::string const &value_snc_type)
    {
        auto dims(get_or_add_dims(ncio, {vname + ".size"}, {count.size()}));

        auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        get_or_put_att(info_v, ncio.rw, "nnz", &_nnz, 1);

        ncio_vector(ncio, count, true, vname+".count", snc_type, dims);
        ncio_vector(ncio, value, true, vname+".value", snc_type, dims);
    }


// ==================================================================


template<class CountT, class IndexT, class ValueT, int RANK>
RLSparseArray_Generator<CountT,IndexT,ValueT,RANK>::
    RLSparseArray_Generator(RLSparseArray const &arr) :
        value_ii(values_generator())
    {
        for (int i=0; i<RANK; ++i) indices_ii.push_back(indices[i].generator());
    }


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
        for (int i=0; i<RANK; ++i) {
            indices[i].ncio(ncio,
                strprintf("%s.indices_%d", vname.c_str(), i),
                count_snc_type, index_snc_type);
        }
        values.ncio(ncio, vname+".values",
            count_snc_type, value_snc_type);
    }


template<class CountT, class IndexT, class ValueT, int RANK>
typename RLSparseArray<CountT,IndexT,ValueT,RANK>::accum_type
RLSparseArray<CountT,IndexT,ValueT,RANK>::
    accum(ValueEqualT &&value_eq)
    {
        std::vector<vaccum_IndexT> vaccum_indices;
        for (int i=0; i<RANK; ++i) {
            vaccum_indices.push_back(indices[i].vaccum(RLAlgo::DIFFS));
        }

        return accum_type(std::move(vaccum_indices),
            values.vaccum(RLAlgo::PLAIN, std::move(value_eq)),
            _nnz);
    }

// ==================================================================








// ---------------------------------------------------------------------
} // namespace ibmisc
#endif
