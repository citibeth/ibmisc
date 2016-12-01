/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IBMISC_COMPACT_MAP
#define IBMISC_COMPACT_MAP

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/iter.hpp>


namespace ibmisc {


// -----------------------------------------
template<class KeyT>
class IndexSet;

// -----------------------------------------

/** Maps keys to integers, according to order of insertion. */
template<class KeyT>        // Our metadata structure
class IndexSet
{
    std::map<KeyT, size_t> _key_to_ix;
    std::vector<KeyT> _ix_to_key;

public:
    IndexSet() {}

    void clear() {
        _key_to_ix.clear();
        _ix_to_key.clear();
    }

    /** Initialize from initializer list */
    IndexSet(std::vector<KeyT> &&keys) : _ix_to_key(std::move(keys))
    {
        for (size_t i=0; i<_ix_to_key.size(); ++i)
            _key_to_ix.insert(std::make_pair(_ix_to_key[i], i));
    }

    std::vector<KeyT> const &keys() const { return _ix_to_key; }

    // --------------------------------------------------------
    typedef EnumRandomIter<typename std::vector<KeyT>::iterator> iterator;
    iterator begin()
        { return enum_random_iter(_ix_to_key.begin(), _ix_to_key.begin()); }
    iterator end()
        { return enum_random_iter(_ix_to_key.begin(), _ix_to_key.end()); }

    typedef EnumRandomIter<typename std::vector<KeyT>::const_iterator> const_iterator;
    const_iterator begin() const
        { return enum_random_iter(_ix_to_key.begin(), _ix_to_key.begin()); }
    const_iterator end() const
        { return enum_random_iter(_ix_to_key.begin(), _ix_to_key.end()); }
    const_iterator cbegin() const
        { return enum_random_iter(_ix_to_key.cbegin(), _ix_to_key.cbegin()); }
    const_iterator cend() const
        { return enum_random_iter(_ix_to_key.cbegin(), _ix_to_key.cend()); }
    // --------------------------------------------------------


    size_t size() const { return _ix_to_key.size(); }

    size_t insert(KeyT const &key)
    {
        if (_key_to_ix.find(key) != _key_to_ix.end()) {
            std::ostringstream buf;
            buf << "Duplicate key detected trying to insert " << key;
            (*ibmisc_error)(-1, "%s", buf.str().c_str());
        }
        _key_to_ix.insert(std::make_pair(key, _key_to_ix.size()));
        size_t ix = _ix_to_key.size();
        _ix_to_key.push_back(key);
        return ix;
    }

    bool contains(KeyT const &key) const
    {
        auto ii(_key_to_ix.find(key));
        return (ii != _key_to_ix.end());
    }


    size_t at(KeyT const &key) const
    {
        auto ii(_key_to_ix.find(key));
        if (ii == _key_to_ix.end()) {
            std::stringstream buf;
            buf << "Cannot find key: " << key;
            (*ibmisc_error)(-1, "%s", buf.str().c_str());
        }
        return ii->second;
    }

    KeyT const &operator[](size_t ix) const
    {
        if (ix >= _ix_to_key.size()) {
            (*ibmisc_error)(-1,
                "Index out of range [0,%ld): %ld", _ix_to_key.size(), ix);
        }
        return _ix_to_key[ix];
    }
};


template<class KeyT>
std::ostream &operator<<(std::ostream &out, ibmisc::IndexSet<KeyT> const &con)
{
    out << "IndexSet([" << std::endl;
    for (size_t i=0; i < con.size(); ++i)
        out << con.at(i) << ", ";
    out << "])";
    return out;
}

}   // Namespace



#endif  // Guard

