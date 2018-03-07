#ifndef IBMISC_UNORDERED_MAP_HPP
#define IBMISC_UNORDERED_MAP_HPP

#include <unordered_map>

namespace ibmisc {

std::unordered_map<class KeyT, class ValueT>
std::unorderd_map<KeyT, ValueT> make_unorderd_map(
    std::vector<std::pair<KeyT,ValueT>> const &inits);

std::unordered_map<class KeyT, class ValueT>
std::unorderd_map<KeyT, ValueT> make_unorderd_map(
    std::vector<std::pair<KeyT,ValueT>> const &inits)
{
    std::unordered_map<KeyT,ValueT> ret;
    for (auto const &pair : inits) {
        ret.insert(pair);
    }
    return ret;
}

}

#endif // IBMISC_UNORDERED_MAP_HPP
