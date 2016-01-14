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
	ValueT operator*() const
		{ return derived()->operator[](0); }

	ValueT *operator->() const
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
template<class ValueT, class SubIterT>
class DerefRandomAccessIter : public random_access_iterator<ValueT, DerefRandomAccessIter<ValueT, SubIterT>>
{
	SubIterT ii;
public:
	DerefRandomAccessIter(SubIterT &&_ii) : ii(_ii) {}

	ValueT &operator[](int n) const
		{ return *ii[n]; }
	DerefRandomAccessIter &operator+=(int n)
		{ii += n; return *this;}
	DerefRandomAccessIter operator+(int n) const
		{ return DerefRandomAccessIter(ii + n); }


	bool operator==(const DerefRandomAccessIter& rhs) const
		{return ii == rhs.ii;}
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
template<class KeyT, class ValueT, class SubIterT>
class SecondIter : public forward_iterator<ValueT, SecondIter<KeyT, ValueT, SubIterT>>
{
	SubIterT ii;
public:
	SecondIter(SubIterT &&_ii) : ii(_ii) {}

	KeyT const &key() const
		{ return ii->first; }

	ValueT &operator*() const
		{ return ii->second; }
	SecondIter &operator++()
		{ ++ii; return *this; }

	bool operator==(const SecondIter& rhs) const
		{return ii == rhs.ii;}
};
// ------------------------------------------------------------



}
#endif	// Guard
