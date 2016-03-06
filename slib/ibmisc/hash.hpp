#pragma once

#include <functional>
#ifdef USE_BOOST
#include <boost/functional/hash.hpp>
#endif

#if 0
namespace ibmisc {

/** Utility class to hash std::pair, which was left out of std::hash.
Used to construct std::unordered_maps where std::pair is used as the
key.  Might have been made obsolete by an implementation of
std::hash<> for std::pair.
@see gridutil.hpp,std::hash< pair< S, T > > */
template<class A, class B>
struct HashPair {
    size_t operator()(std::pair<A, B> const &x) const throw()
        { return std::hash<A>()(x.first)*31 + std::hash<B>()(x.second); }
};

}
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
