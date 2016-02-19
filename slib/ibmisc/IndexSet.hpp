#ifndef IBMISC_COMPACT_MAP
#define IBMISC_COMPACT_MAP

#include <ibmisc/ibmisc.hpp>
#include <vector>
#include <map>

namespace ibmisc {


// -----------------------------------------
template<class KeyT>
class IndexSet;

template<class KeyT>
std::ostream &operator<<(std::ostream &out, IndexSet<KeyT> const &con);
// -----------------------------------------

/** Maps keys to integers, according to order of insertion. */
template<class KeyT>		// Our metadata structure
class IndexSet
{
	std::map<KeyT, size_t> _key_to_ix;
	std::vector<KeyT> _ix_to_key;

public:
	typedef typename std::vector<KeyT>::iterator iterator;
	iterator begin() { return _ix_to_key.begin(); }
	iterator end() { return _ix_to_key.end(); }

	typedef typename std::vector<KeyT>::const_iterator const_iterator;
	const_iterator begin() const { return _ix_to_key.begin(); }
	const_iterator end() const { return _ix_to_key.end(); }
	const_iterator cbegin() const { return _ix_to_key.cbegin(); }
	const_iterator cend() const { return _ix_to_key.cend(); }


	size_t size() { return _ix_to_key.size(); }

	void insert(KeyT const &key)
	{
		if (_key_to_ix.find(key) != _key_to_ix.end()) {
			std::stringstream buf;
			buf << "Duplicate key detected trying to insert " << key;
			(*ibmisc_error)(-1, "%s", buf.str().c_str());
		}
		_key_to_ix.insert(std::make_pair(key, _key_to_ix.size()));
		_ix_to_key.push_back(key);
	}

	bool contains(KeyT const &key)
	{
		auto ii(_key_to_ix.find(key));
		return (ii != _key_to_ix.end());
	}

	bool contains(size_t ix)
	{
		return (ix < _ix_to_key.size());
	}


	size_t at(KeyT const &key)
	{
		auto ii(_key_to_ix.find(key));
		if (ii == _key_to_ix.end()) {
			std::stringstream buf;
			buf << "Cannot find key: " << key;
			(*ibmisc_error)(-1, "%s", buf.str().c_str());
		}
		return ii->second;
	}

	KeyT const &operator[](size_t ix)
	{
		if (ix >= _ix_to_key.size()) {
			(*ibmisc_error)(-1,
				"Index out of range [0,%ld): %ld", _ix_to_key.size(), ix);
		}
		return _ix_to_key[ix];
	}

	// http://www.cplusplus.com/forum/general/45776/
	friend std::ostream &operator<< <>(std::ostream &out, IndexSet<KeyT> const &con);

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


}	// Namespace


#endif	// Guard
