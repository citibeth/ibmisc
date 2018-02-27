#ifndef SPSPARSE_SPARSEARRAY_HPP
#define SPSPARSE_SPARSEARRAY_HPP

namespace spsparse {
namespace accum {

/** General Accumulator for spare arrays; aggreagates vaccumulators for each index and value. */
template<class IndexVAccumT, class ValueVAccumT, int RANK>
class SparseArrayAgg {

    std::vector<IndexVAccumT> vaccum_index;
    ValueVAccumT vaccum_value;

    SparseArrayAgg(
        std::vector<IndexVAccumT> &&_vaccum_index,
        ValueVAccumT &&_vaccum_value)
    : vaccum_index(std::move(_vaccum_index)), vaccum_value(std::move(_vaccum_value)) {}

    inline void add(std::array<index_type,RANK> const &index, val_type const &val)
    {
        for (int i=0; i<RANK; ++i) vaccum_index[i].add(index[i]);
        vaccum_value.add(val);
    }
}

}}    // namespace accum
#endif    // guard
