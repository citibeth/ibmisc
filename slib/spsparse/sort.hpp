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

#pragma once

#include <vector>
#include <algorithm>
#include <blitz/array.h>
#include <spsparse/spsparse.hpp>

/** This is held over from the "old" giss library.  It is used in
IceBin. The alternative, to create a Blitz-based CooVector, is more
work. */

namespace spsparse {

namespace _sorted_perm_private {

template<class KeyT>
struct CmpIndex1 {
    blitz::Array<int, 1> keys1;
public:
    void init(blitz::Array<int, 1> const &_keys1) {
        keys1.reference(_keys1);
    }
    bool operator()(int i, int j)
    {
        return (keys1(i) < keys1(j));
    }
};
}
    

/** Creates a permutation that, if applied, will result in a sorted
version of the array.  The user is responsible for applying it as
needed, later.
@param ikeys The array (of keys) to sort
@return The permutation, stored in a std::vector.  sorted_ikeys =
ikeys(perm[i])) will produce the sorted vector.
*/
template<class KeyT>
std::vector<int> sorted_perm(
    blitz::Array<KeyT,1> const &ikeys)
{

    _sorted_perm_private::CmpIndex1<KeyT> cmp;
    cmp.init(ikeys);

    // Generate a permuatation
    int n = ikeys.size();
    std::vector<int> perm; perm.reserve(n);
    for (int i=0; i<n; ++i) perm.push_back(i);
    std::sort(perm.begin(), perm.end(), cmp);

    return perm;
}

//enum class DuplicatePolicy {REPLACE, ADD};

/** Combination of consolidate duplicate indices, and apply permutation.
Consolidation means removing duplicate entries from a sparse vector representation.
@param The original array of indices.
@param A permutation that will sort indices (output of sorted_perm()).
@param ivalues Array of values corresponding to the indices array.
@param ovalues OUTPUT: Array of consolidated (must already be allocated).
@param duplicate_policy Upon finding duplicate indices, should we sum
or replace in the values array?  When consolidating values of a sparse
vector, usually use ADD; when consolidating indices, use REPLACE*/
template<class IndexT, class ValT>
int consolidate_by_perm(
    blitz::Array<IndexT,1> const &indices,
    std::vector<int> const &sorted_perm,
    blitz::Array<ValT,1> const &ivalues,
    blitz::Array<ValT,1> &ovalues,
    DuplicatePolicy duplicate_policy)
{
    // Location to put this in the final array
    std::vector<int> dest; dest.reserve(sorted_perm.size());

    // Scan through, overwriting our array
    // New array will never be bigger than original
    int j=0;        // Last-written item in output array
    IndexT last_index = indices(sorted_perm[0]);    // Last index we saw in input array
    for (unsigned int i=1; i<sorted_perm.size(); ++i) {
        IndexT const &index = indices(sorted_perm[i]);
        ValT const &val = ivalues(sorted_perm[i]);

        if (index == last_index) {
            if (duplicate_policy == DuplicatePolicy::ADD)
                ovalues(j) += val;
            else
                ovalues(j) = val;
        } else {
            ++j;
            ovalues(j) = val;
            last_index = index;
        }
    }
    int n = j+1;    // Size of output array

    return n;
}

}
