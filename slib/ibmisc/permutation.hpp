#include <algorithm>

namespace ibmisc {

template<class RandomIt, class IndexT>
void sorted_permutation(RandomIt in0, RandomIt in1, std::vector<IndexT> &perm)
{
    size_t const N = in1 - in0;

    // Construct permutation
    IndexT i = 0;
    for (RandomIt ii(in0); ii != in1; ++ii, ++i)
        perm.push_back(i);

    // Sort it
    std::sort(perm.begin(), perm.end(),
        [&](IndexT a, IndexT b)
            { return in0[a] < in0[b]; }
    );
}

}    // namespace
