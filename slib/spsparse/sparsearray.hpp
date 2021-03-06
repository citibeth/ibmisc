#ifndef SPSPARSE_SPARSEARRAY_HPP
#define SPSPARSE_SPARSEARRAY_HPP

#include <array>

namespace spsparse {
namespace accum {

/** General Accumulator for spare arrays; aggreagates vaccumulators for each index and value. */
template<class IndexVAccumT, class ValueVAccumT, int RANK>
class SparseArrayAgg {

    std::vector<IndexVAccumT> vaccum_index;
    ValueVAccumT vaccum_value;
    int &nnz;    // Used to count number of elements added

public:
    typedef typename IndexVAccumT::val_type index_type;
    typedef typename ValueVAccumT::val_type value_type;

    SparseArrayAgg(
        std::vector<IndexVAccumT> &&_vaccum_index,
        ValueVAccumT &&_vaccum_value,
        int &_nnz)
    : vaccum_index(std::move(_vaccum_index)),
      vaccum_value(std::move(_vaccum_value)),
      nnz(_nnz) {}

    inline void add(std::array<index_type, RANK> const &index, value_type const &val)
    {
        for (int i=0; i<RANK; ++i) vaccum_index[i].add(index[i]);
        vaccum_value.add(val);
        ++nnz;
    }
};

}}    // namespace accum
#endif    // guard
