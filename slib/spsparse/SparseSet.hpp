#pragma once

#include <unordered_map>
#include <algorithm>
#include <spsparse/spsparse.hpp>

namespace spsparse {

/** Translates between a sparse set (say, the set of indices used in a
SparseVector) and a dense set numbered [0...n) */
template<class SparseT, class DenseT>
class SparseSet {
	SparseT _sparse_extent;
	std::unordered_map<SparseT, DenseT> _s2d;
	std::vector<SparseT> _d2s;

public:
	SparseSet() : _sparse_extent(-1) {}

	bool in_sparse(SparseT sparse_ix)
		{ return _s2d.find(sparse_ix) != _s2d.end(); }

	SparseT sparse_extent() const
		{ return _sparse_extent; }

	void set_sparse_extent(SparseT extent)
		{ _sparse_extent = extent; }

	DenseT dense_extent() const
		{ return _d2s.size(); }

	void add(SparseT sparse_index)
	{
		if (_s2d.find(sparse_index) == _s2d.end()) {
			_s2d.insert(std::pair<SparseT,DenseT>(sparse_index, dense_extent()));
			_d2s.push_back(sparse_index);
		}
	}

	template<class IterT>
	void add(IterT sparsei, IterT const &sparse_end)
	{
		for (; sparsei != sparse_end; ++sparsei) add(*sparsei);
	}


	template<class IterT>
	void add_sorted(IterT sparsei, IterT const &sparse_end)
	{
		std::vector<SparseT> sorted;
		for (; sparsei != sparse_end; ++sparsei) sorted.push_back(*sparsei);
		std::sort(sorted.begin(), sorted.end());

		add(sorted.begin(), sorted.end());
	}

	DenseT to_dense(SparseT const &sval) const
		{ return _s2d.at(sval); }

	SparseT to_sparse(DenseT const &dval) const
	{
		if (dval < 0 || dval >= _d2s.size()) {
			(*spsparse_error)(-1,
				"Value %ld is out of range (0, %ld)", dval, _d2s.size());
		}
		return _d2s[dval];
	}
};

}	// Namespace
