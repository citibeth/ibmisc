#ifndef IBMISC_ITERATOR_HPP
#define IBMISC_ITERATOR_HPP

namespace ibmisc {

/**An interator that derferences its result one more time than the
wrapped iterator.  Implemented as a wrapper.  Useful to iterate
through colletions of unique_ptr as if they're plain values.  For
example:
@code
 MyIterator ii = ...;
 DerefIterator<MyIterator> dii(ii);
 **ii == *dii
@endcode
@param IterT The type of iterator to wrap.
@see RerefIterator
*/
template<class IterT>
class DerefIter : public IterT {
public:
	/** Construct by wrapping an existing iterator */
	DerefIterator(IterT const &ii) : IterT(ii) {}

	auto operator*() -> decltype(*(this->IterT::operator*()))
		{ return *(this->IterT::operator*()); }
	auto operator->() -> decltype(&*(this->IterT::operator*()))
		{ return &*(this->IterT::operator*()); }

};

template<class IterT>
inline DerefIter<IterT> deref_iter(IterT &&ii)
	{ return DerefIter<IterT>(ii); }


}
#endif	// Guard
