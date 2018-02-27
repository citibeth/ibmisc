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

    void ncio(NcIO &ncio, std::string const &vname, std::string const &count_snc_type, std::string const &value_snc_type)
    {
        auto dims(get_or_add_dims(ncio, {vname + ".size"}, {count.size()}));
        ncio_vector(ncio, count, true, vname+".count", snc_type, dims);
        ncio_vector(ncio, value, true, vname+".value", snc_type, dims);
    }

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

/** NOTE: Singular used for template parameters, plural used for
    vectors holding plural things (or accumulators). */
template<class CountT, class IndexT, class ValueT, int RANK>
struct RLSparseArray {

    typedef RLVector<CountT,IndexT> indices_rlvector_type;
    typedef RLVector<CountT,ValueT> values_rlvector_type;

    typedef  accum_type;

    std::array<indices_rlvector_type,RANK> indices;

    void ncio(NcIO &ncio, std::string const &vname,
        std::string const &count_snc_type,
        std::string const &index_snc_type,
        std::string const &value_snc_type)
    {
        for (int i=0; i<RANK; ++i) {
            indices[i].ncio(ncio,
                strprintf("%s.indices_%d", vname.c_str(), i),
                count_snc_type, index_snc_type);
        }
        values.ncio(ncio, vname+".indices",
            count_snc_type, value_snc_type);
    }

    /** Creates an encoder */
    template<class ValueEqualT>
    accum::SparseArrayAgg<
        typename indices_rlvector_type::vaccum_type,
        typename values_rlvector_type::vaccum_type,
        RANK>
    accum(ValueEqualT &&value_eq) const
    {
        std::vector<vaccum_IndexT> vaccum_indices;
        for (int i=0; i<RANK; ++i) {
            vaccum_indices.push_back(indices[i].vaccum(RLAlgo::DIFFS));
        }

        return AccumT(std::move(vaccum_indices),
            values.vaccum(RLAlgo::PLAIN, std::move(value_eq)));
    }

    class Generator {
        std::vector<typename indices_rlvector_type::generator_type> indices_ii;
        typename values_rlvector_type::generator_type value_ii;

        Generator(RLSparseArray const &arr) :
            value_ii(values_generator())
        {
            for (int i=0; i<RANK; ++i) indices_ii.push_back(indices[i].generator());
        }


        bool operator++()
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

        IndexT index(int ix) const
            { return *indices_ii[ix]; }
        std::array<IndexT,RANK> index() const
        {
            std::array<IndexT,RANK> ret;
            for (int i=0; i<RANK; ++i) ret[i] = index(i);
            return ret;
        }

        double value() const
            { return *values_ii; }

    };

    Generator generator() const
        { return Generator(*this); }


    template<class AccumT>
    void spcopy(AccumT &&ret)
    {
        for (auto gen(Generator); ++gen;) {
            ret.add(gen.index(), gen.value());
        }
    }


};


TODO: This goes into IceBin

template<class CountT, class int, class ValueT>
void weighted_sparse_apply(
    blitz::Array<double,1> &B,
    blitz::Array<double,1> const &wM,            // ==0 in M's nullspace
    RLSparseArray<int,int,double,2> const &M,    // BvA
    blitz::Array<double,1> const &Mw,            // ==0 in M's nullspace
    blitz::Array<double,1> const &A,
    bool clear,
    double fill_value,
    bool check_for_fill,
    bool force_conservation)
{
    // Clear the output
    if (clear) B = fill_value;

    // Multiply
    for (auto gen(M.generator()); ++gen; ) {
        auto const iB(gen.index(0));
        auto const iA(gen.index(1));
        double const val = A(iA) * gen.value();
        if (check_for_fill) {
            if (B(iB) == fill_value) B(iB) = val;
            else B(iB) += val;
        } else {
            B(iB) += val;
        }
    }

    // Adjust for conservation
    if (force_conservation && !M_is_conservative) {
        // Determine original mass
        double mass_A = 0;
        for (int i=0; i<A.extent(0); ++i) {
            if (Mw(i) != 0) mass_A += Mw(i) * A(i);
        }

        // Determine final mas
        double mass_B = 0;
        for (int i=0; i<B.extent(0); ++i) {
            if (wM(i) != 0) mass_B += Mw(i) * B(i);
        }

        // Determine adjustment factor for B
        double const factor = mass_A / mass_B;

        // Multiply by the adjustment factor
        for (int i=0; i<B.extent(0); ++i) {
            if (wM(i) != 0) B(i) *= factor;
        }
    }

}






// ---------------------------------------------------------------------
} // namespace ibmisc
#endif
