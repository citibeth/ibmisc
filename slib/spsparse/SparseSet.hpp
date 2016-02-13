#pragma once

#include <unordered_map>
#include <algorithm>
#include <spsparse/spsparse.hpp>

namespace spsparse {

/** Translates between a sparse set (say, the set of indices used in a
SparseVector) and a dense set numbered [0...n) */
template<class SparseT, class DenseT>
class SparseSet {
	std::unordered_map<SparseT, DenseT> _s2d;
	std::vector<SparseT> _d2s;

public:
	size_t size() { return _d2s.size(); }

	template<class IterT>
	void add(IterT sparsei, IterT const &sparse_end)
	{
		for (; sparsei != sparse_end; ++sparsei) {
			if (_s2d.find(*sparsei) == _s2d.end()) {		// Add a new index
				_s2d.insert(std::pair<SparseT,DenseT>(*sparsei, size()));
				_d2s.push_back(*sparsei);
			}
		}
	}


	template<class IterT>
	void add_sorted(IterT sparsei, IterT const &sparse_end)
	{
		std::vector<SparseT> sorted;
		for (; sparsei != sparse_end; ++sparsei) sorted.push_back(*sparsei);
		std::sort(sorted.begin(), sorted.end());

		add(sorted.begin(), sorted.end());
	}

	DenseT to_dense(SparseT const &sval)
		{ return _s2d.at(sval); }

	SparseT to_sparse(DenseT const &dval)
	{
		if (dval < 0 || dval >= _d2s.size()) {
			(*spsparse_error)(-1,
				"Value %ld is out of range (0, %ld)", dval, _d2s.size());
		}
		return _d2s[dval];
	}
};

}	// Namespace
