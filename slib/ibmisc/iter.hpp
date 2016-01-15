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



}
#endif	// Guard
