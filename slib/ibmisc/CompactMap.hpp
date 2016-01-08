#ifndef IBMISC_COMPACT_MAP
#define IBMISC_COMPACT_MAP

#include <ibmisc/ibmisc.hpp>
#include <vector>
#include <map>

namespace ibmisc {


// -----------------------------------------
template<class KeyT, class ValueT>
class CompactMap;

template<class KeyT, class ValueT>
std::ostream &operator<<(std::ostream &out, CompactMap<KeyT,ValueT> const &con);
// -----------------------------------------

/** A two-way mapping between strings and ints (densely spaced), that
can be changed at runtime.  This serves as a virtual base class,
allowing algorithms to be written in terms of DynamicEnum (eg CouplingContract).

The string "unit" is special.  By convention, it is the LAST item in the enum.
Size can be queried with or without the "unit" string included.


ValueT must provide KeyT ValueT::name()
(This could be made more general with templated accessors.  Ugh)
*/
template<class KeyT, class ValueT>		// Our metadata structure
class CompactMap
{
	std::vector<ValueT> _ix_to_value;
	std::map<KeyT, size_t> _key_to_ix;

public:
	CompactMap() {}

	// ----------------------------------------------
	typedef typename std::vector<ValueT>::const_iterator iterator;
	typedef typename std::vector<ValueT>::const_iterator const_iterator;

	iterator begin()
		{ return _ix_to_value.begin(); }
	iterator end()
		{ return _ix_to_value.end(); }

	const_iterator begin() const
		{ return _ix_to_value.begin(); }
	const_iterator end() const
		{ return _ix_to_value.end(); }
	// ----------------------------------------------

	size_t insert(ValueT &&cf);
	size_t insert(ValueT const &cf);
	size_t size() const { return _ix_to_value.size(); }
	int index(std::string const &name, bool raise_error=true) const;

	ValueT const &operator[](int ix) const
		{ return _ix_to_value[ix]; }

	/** Our behavior here is closest to STL's std::map<>::at(). */
	ValueT const &operator[](std::string const &name) const
		{ return (*this)[this->index(name)]; }



	// http://www.cplusplus.com/forum/general/45776/
	friend std::ostream &operator<< <>(std::ostream &out, CompactMap<KeyT,ValueT> const &con);

};


template<class KeyT, class ValueT>
std::ostream &operator<<(std::ostream &out, ibmisc::CompactMap<KeyT,ValueT> const &con)
{
	out << "CompactMap(" << std::endl;

	for (size_t i=0; i < con.size(); ++i)
		out << "    " << i << ": " << con[i].name() << std::endl;

	out << ")";
	return out;
}



// ------------------------------------------------------------

template<class KeyT, class ValueT>
	size_t CompactMap<KeyT,ValueT>::insert(ValueT &&cf)
	{
		size_t index = _ix_to_value.size();
		_key_to_ix.insert(std::make_pair(cf.name(), index));
		_ix_to_value.push_back(std::move(cf));
		return index;
	}

template<class KeyT, class ValueT>
	size_t CompactMap<KeyT,ValueT>::insert(ValueT const &cf)
	{
		size_t index = _ix_to_value.size();
		_key_to_ix.insert(std::make_pair(cf.name(), index));
		_ix_to_value.push_back(cf);
		return index;
	}

template<class KeyT, class ValueT>
	int CompactMap<KeyT,ValueT>::index(std::string const &name, bool raise_error) const
	{
		auto ii = _key_to_ix.find(name);
		if (ii == _key_to_ix.end()) {
			if (raise_error) {
				(*ibmisc_error)(-1,
					"CompactMap::index(): key not found");
			} else return -1;
		}
		return ii->second;
	}



}	// Namespace


#endif	// Guard
