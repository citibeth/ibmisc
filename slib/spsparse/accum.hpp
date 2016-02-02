#ifndef SPSPARSE_ACCUM_HPP
#define SPSPARSE_ACCUM_HPP

#include <cstddef>
#include <array>
#include <vector>
#include <ibmisc/blitz.hpp>
#include <spsparse/spsparse.hpp>

namespace spsparse {

/** @defgroup accum accum.hpp
@brief Accumulators to use with algorithms.

An accumulator is used by an algorithm to accept output via an add()
method.  This allows each algorithm to support a variety of output
strategies: in-place, out-of-place, conversion to dense arrays, sum
for inner product, etc.

TODO: Figure out an implement good semantics for set_shape().  When is
it called / not called?

@{
*/

// -----------------------------------------------------------
/** @brief For in-place operations.

Overwrites a VectorCooArray::iterator.  Good for
in-place operations that we KNOW won't change the size of the
VectorCooArray (eg: use with spsparse::transpose()).

@warning Does not bounds-check for the END of the iterator.

Usage Example:
@code
VectorCooArray<int,double,2> A;
OverWriteAccum<decltype(A)::iterator> overwrite(A.begin());
transpose(overwrite, A, {1,0});
@endcode

*/
template<class IterT>
class OverwriteAccum
{
	SPSPARSE_LOCAL_TYPES(IterT);

	IterT ii;
public:
	OverwriteAccum(IterT &&_ii) : ii(std::move(_ii)) {}

	void add(indices_type const &index, typename IterT::val_type const &val) {
		ii.set_index(index);
		ii.val() = val;
		++ii;
	}
};
// -----------------------------------------------------------
/** @brief For transpose or project.

Permutes/selects dimensions, storing in a sub-accumulator.  The
sub-accumulator does NOT have to have the same rank; this may
therefore be used to transpose or remove dimensions.

Usage Example:
@code
VectorCooArray<int,double,2> A,B;
PermuteAccum<A::rank, decltype(B)> p(B, {1,0});
copy(p, A);
@endcode

*/
template<int IN_RANK, class AccumulatorT>
class PermuteAccum
{
public:
	static const int rank = IN_RANK;
	static const int out_rank = AccumulatorT::rank;

private:
	AccumulatorT sub;
	std::vector<int> perm;
	std::array<int, out_rank> out_idx;

public:
	PermuteAccum(AccumulatorT &&_sub, std::vector<int> const &_perm)
		: sub(std::move(_sub)), perm(_perm) {}

	void add(std::array<int, IN_RANK> const &index, typename AccumulatorT::val_type const &val) {
		for (int i=0; i<IN_RANK; ++i) out_idx[i] = index[perm[i]];
		sub.add(out_idx, val);
	}
};
// -----------------------------------------------------------
/** @brief For output/conversion to dense arrays.

Outputs to a blitz::Array.  Note that blitz:Array is quite flexible,
and can be used to access almost any existing fixed-size data
structure.

Usage Example:
@code
VectorCooArray<int,double,2> A;
blitz::Array<double,2> B(to_tiny(A.shape));
DenseAccum<decltype(A)> Baccum(B);
copy(Baccum, A);
@endcode

*/
template<class VectorCooArrayT>
struct DenseAccum
{
	SPSPARSE_LOCAL_TYPES(VectorCooArrayT);

private:
	DuplicatePolicy duplicate_policy;
	blitz_type dense;
	blitz::TinyVector<int, rank> bidx;

public:
	DenseAccum(blitz_type &_dense, DuplicatePolicy _duplicate_policy=DuplicatePolicy::ADD) : dense(_dense), duplicate_policy(_duplicate_policy) {}

	void add(indices_type const &index, val_type const &val)
	{
		ibmisc::to_tiny<int, index_type, rank>(bidx, index);
		val_type &oval(dense(bidx));

		switch(duplicate_policy) {
			case DuplicatePolicy::LEAVE_ALONE :
				if (!std::isnan(oval)) oval = val;
			break;
			case DuplicatePolicy::ADD :
				oval += val;
			break;
			case DuplicatePolicy::REPLACE :
				oval = val;
			break;
		}
	}
};
// -----------------------------------------------------------


// -------------------------------------------------------
/** @brief For inner product.

Accumulates into a single scalar, ignoring index information.  Used to
implement inner products, or just sum over an array.

Usage Example:
@code
VectorCooArray<int,double,2> A;
ScalarAccumulator<decltype(A)> s;
copy(s, A);
printf("Sum = %g\n", s.val);
@endcode
*/
template<class VectorCooArrayT>
struct ScalarAccumulator {
	SPSPARSE_LOCAL_TYPES(VectorCooArrayT);

	val_type val;
	ScalarAccumulator() : val(0) {}

	void add(const std::array<index_type, rank> &index, val_type const &_val)
		{ this->val += _val; }
};




/** @} */


}	// Namespace

#endif // Guard
