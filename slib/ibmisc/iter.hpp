/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IBMISC_ITERATOR_HPP
#define IBMISC_ITERATOR_HPP

#include <iterator>

namespace ibmisc {


/** Base class with useful stuff in it.

http://www.cplusplus.com/reference/iterator/RandomAccessIterator/
https://groups.google.com/a/isocpp.org/forum/#!msg/std-proposals/tuYkpJqCKDA/r_1c5V97MckJ

copy and move constructors
X a;
X b(a);
b = a;
a == b
operator[]
operator+
operator+=

https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

 */
template<class ValueT, class DerivedT>
class random_access_iterator : public std::iterator<std::random_access_iterator_tag, ValueT>
{
    DerivedT *derived()
        { return static_cast<DerivedT *>(this); }

    DerivedT const *derived() const
        { return static_cast<DerivedT const *>(this); }

public:
    ValueT &operator*()
        { return derived()->operator[](0); }

    ValueT *operator->()
        { return &operator*(); }

    DerivedT& operator++()
        { return derived()->operator+=(1); }

#if 1
    DerivedT operator++(int)
    {
        DerivedT clone(derived());
        derived()->operator+=(1);
        return clone;
    }
#endif

    DerivedT& operator--()
        { return derived()->operator+=(-1); }
#if 1
    DerivedT operator--(int)
    {
        DerivedT clone(derived());
        derived()->operator+=(-1);
        return clone;
    }
#endif

    DerivedT &operator-=(int n)
        { return derived()->operator+=(-n); }

    bool operator!=(const DerivedT& rhs) const
        {return !derived()->operator==(rhs); }
};

// ------------------------------------------------------------
// ------------------------------------------------------------
// http://www.cplusplus.com/reference/iterator/iterator/
template<class ValueT, class WrappedIterT>
class DerefRandomAccessIter : public random_access_iterator<ValueT, DerefRandomAccessIter<ValueT, WrappedIterT>>
{
public:
    WrappedIterT wrapped;

    DerefRandomAccessIter(WrappedIterT &&ii) : wrapped(ii) {}

    ValueT &operator[](int n) const
        { return *wrapped[n]; }
    DerefRandomAccessIter &operator+=(int n)
        {wrapped += n; return *this;}
    DerefRandomAccessIter operator+(int n) const
        { return DerefRandomAccessIter(wrapped + n); }


    bool operator==(const DerefRandomAccessIter& rhs) const
        {return wrapped == rhs.wrapped;}
};
// ------------------------------------------------------------
template<class ValueT, class DerivedT>
class forward_iterator : public std::iterator<std::forward_iterator_tag, ValueT>
{
    DerivedT *derived()
        { return static_cast<DerivedT *>(this); }
    DerivedT const *derived() const
        { return static_cast<DerivedT const *>(this); }
public:
    ValueT *operator->() const
        { return &derived()->operator*(); }

    DerivedT operator++(int)
    {
        DerivedT clone(derived());
        derived()->operator+=(1);
        return clone;
    }

    bool operator!=(const DerivedT& rhs) const
        {return !derived()->operator==(rhs); }
};
// -------------------------------------------------------
template<class KeyT, class ValueT, class WrappedIterT>
class SecondIter : public forward_iterator<ValueT, SecondIter<KeyT, ValueT, WrappedIterT>>
{
public:
    WrappedIterT wrapped;

    SecondIter(WrappedIterT &&ii) : wrapped(ii) {}

    KeyT const &key() const
        { return wrapped->first; }

    ValueT &operator*() const
        { return wrapped->second; }
    SecondIter &operator++()
        { ++wrapped; return *this; }

    bool operator==(const SecondIter& rhs) const
        {return wrapped == rhs.wrapped;}
};
// ------------------------------------------------------------
template<class KeyT, class ValueT, class WrappedIterT>
class DerefSecondIter : public forward_iterator<ValueT, DerefSecondIter<KeyT, ValueT, WrappedIterT>>
{
public:
    WrappedIterT wrapped;
    DerefSecondIter(WrappedIterT &&ii) : wrapped(ii) {}

    KeyT const &key() const
        { return wrapped->first; }

    ValueT &operator*() const
        { return *wrapped->second; }

    DerefSecondIter &operator++()
        { ++wrapped; return *this; }

    bool operator==(const DerefSecondIter& rhs) const
        {return wrapped == rhs.wrapped;}
};
// ------------------------------------------------------------
template<class IterT>
class EnumRandomIter
{
public:
    typedef typename IterT::value_type value_type;

protected:
    IterT const begin;
    IterT ii;

public:
    EnumRandomIter(IterT const &_begin, IterT const &_ii) : begin(_begin), ii(_ii) {}


    value_type operator[](int n)
        { return ii[n]; }
    value_type operator*()
        { return *ii; }

    EnumRandomIter &operator+=(int n)
        { ii += n; return *this; }
    EnumRandomIter& operator--()
        { --ii; return *this; }
    EnumRandomIter& operator++()
        { ++ii; return *this; }
    EnumRandomIter &operator-=(int n)
        { ii -= n; return *this; }

    EnumRandomIter operator+(int n) const
        { return EnumRandomIter(begin, ii+n); }

    bool operator==(EnumRandomIter const &rhs) const
        { return ii == rhs.ii; }
    bool operator!=(const EnumRandomIter& rhs) const
        {return ii != rhs.ii; }

#if 0
    bool operator==(IterT const &rhs) const
        { return ii == rhs; }
    bool operator!=(const IterT& rhs) const
        {return ii != rhs; }
#endif

    size_t index()
        { return ii - begin; }
};

/** Convenience constructor */
template<class IterT>
EnumRandomIter<IterT> enum_random_iter(IterT const &_begin, IterT const &_ii)
    { return EnumRandomIter<IterT>(_begin, _ii); }


// ------------------------------------------------------------



}
#endif  // Guard
