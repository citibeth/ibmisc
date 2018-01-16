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

#include <functional>
#ifdef USE_BOOST
#include <boost/functional/hash.hpp>
#endif

/**@file std::hash<std::pair<>> is defined in "namespace std" to fix a shortcoming in C++.

@see: http://stackoverflow.com/questions/7222143/unordered-map-hash-function-c
*/
namespace std
{
  /** Used to hash elements of type std::pair<>.  This should have been
  in the standard C++11
  @see: http://stackoverflow.com/questions/7222143/unordered-map-hash-function-c */
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
#ifdef USE_BOOST
    {
      size_t seed = 0;
      boost::hash_combine(seed, v.first);
      boost::hash_combine(seed, v.second);
      return seed;
    }
#else
    { return std::hash<S>()(v.first)*31 + std::hash<T>()(v.second); }
#endif

  };
}
