#ifndef SPSPARSE_ALGORITHM_HPP
#define SPSPARSE_ALGORITHM_HPP

#include <memory>
#include <spsparse/spsparse.hpp>
#include <spsparse/xiter.hpp>

namespace spsparse {

/** @defgroup algorithm algorithm.hpp
@brief Algorithms on SpSparse arrays

Most operations on SpSpare arrays are written in a type-independent fashion, as algorithms.  Algorithms are all written with input and output parameters, the output being an accumulator.

The first parameter for an algorithm is always the output (an accumulator).  Interestingly, that is always the LAST argument of its template.

@{
*/

// ====================================================


// -----------------------------------------------------------
/** @brief Copy a sparse array
@param ret Accumulator for output.
@param A Input. */
template<class VectorCooArrayT, class AccumulatorT>
void copy(AccumulatorT &ret, VectorCooArrayT const &A);

template<class VectorCooArrayT, class AccumulatorT>
void copy(AccumulatorT &ret, VectorCooArrayT const &A)
{
	std::array<int,VectorCooArrayT::rank> idx;
	for (auto ii=A.begin(); ii != A.end(); ++ii) {
		ret.add(ii.index(), ii.val());
	}
}
// -----------------------------------------------------------
/** @brief Transpose a sparse array (reorder its dimensions)
@param ret Accumulator for output.
@param A Input.
@param perm Permutation on dimensions.  ret.dim[i] == A.dim[perm[i]] */
template<class VectorCooArrayT, class AccumulatorT>
void transpose(AccumulatorT &ret, VectorCooArrayT const &A, std::array<int,VectorCooArrayT::rank> const &perm);

template<class VectorCooArrayT, class AccumulatorT>
void transpose(AccumulatorT &ret, VectorCooArrayT const &A, std::array<int,VectorCooArrayT::rank> const &perm)
{
	std::array<int,VectorCooArrayT::rank> idx;
	for (auto ii=A.begin(); ii != A.end(); ++ii) {
		for (int new_k=0; new_k < VectorCooArrayT::rank; ++new_k) {
			int old_k = perm[new_k];
			idx[new_k] = ii.index(old_k);
		}
		ret.add(idx, ii.val());
	}
}
// -----------------------------------------------------------
/** @brief Determines offset of beginning of each row/col (leading sorted dimension) in an array.
@note The array MUST be sorted properly beforehand, or this will fail.

Code Example
@code
VectorCooMatrix<int, double> A;
A.consolidate({0,1});
std::vector<size_t> row_start = dim_beginnings(A);
A.consolidate({1,0});
std::vector<size_t> col_start = dim_beginnings(A);

See: VectorCooArray::consolidate() */
template<class VectorCooArrayT>
std::vector<size_t> dim_beginnings(VectorCooArrayT const &A);

template<class VectorCooArrayT>
std::vector<size_t> dim_beginnings(VectorCooArrayT const &A)
{
	const int RANK = VectorCooArrayT::rank;

	std::vector<size_t> abegin;

	// Check that we're sorted by SOME dimension.
	if (A.sort_order[0] < 0) {
		(*spsparse_error)(-1, "dim_beginnings() required the VectorCooArray is sorted first.");
	}

	// Get beginning of each row in a (including sentinel at end)
	auto ai(A.begin());
	auto const end(A.end());
	if (ai != end) {		// At least 1 element
		abegin.push_back(ai.offset());			// First item in array is always 0
		int const dim = A.sort_order[0];		// Dimension we're sectioning by
		int last_row = ai.index(dim);

		for (++ai; ; ++ai) {
			if (ai == end) {
				// Add a sential row
				abegin.push_back(ai.offset());
				break;
			}
			if (ai.index(dim) != last_row) {
				// We see a new row!
				abegin.push_back(ai.offset());
				last_row = ai.index(dim);
			}
		}
	}

#if 0
std::cout << "dim_beginnings matrix = " << A << std::endl;
printf("dim_beginnings(size=%ld) = [", A.size());
for (auto ii(abegin.begin()); ii != abegin.end(); ++ii) {
	printf("%d ", *ii);
}
printf("]\n");
#endif

	return abegin;
}
// ----------------------------------------------------------
/** @brief Iterates through a sorted sparse array on a per-row (or column) basis.

Example Code:
@code
VectorCooMatrix<int, double> const A;
VectorCooMatrix<int, double> Arm;
consolidate(Arm, A, {0,1});	// Arm = A row major
auto row_beginnings(dim_beginnings(Arm));
for (DimBeginningsXiter<decltype(A)> ii(&A, 0, 1,
	row_beginnings.begin(), row_beginnings.end());
	!ii.eof(); ++ii)
{
	printf("Row %d contains:", *ii);
	for (auto jj=join_a->sub_xiter(); !jj.eof(); ++jj) {
		printf(" (col=%d : val=%g)", *jj, jj.val());
		}
	printf("\n");
}
@endcode


Convenience methods allow a simplification as follows:
@code
VectorCooMatrix<int, double> A;

A.consolidate({0,1});
for (auto ii(A.dim_beginnings_xiter(); !ii.eof(); ++ii) {
{
	printf("Row %d contains:", *ii);
	for (auto jj=join_a->sub_xiter(); !jj.eof(); ++jj) {
		printf(" (col=%d : val=%g)", *jj, jj.val());
	}
	printf("\n");
}
@endcode

We can also use this to scan through the matrix in row-major order:
@code
VectorCooMatrix<int, double> A;

A.consolidate({1,0});
for (auto ii(A.dim_beginnings_xiter(); !ii.eof(); ++ii) {
{
	printf("Col %d contains:", *ii);
	for (auto jj=join_a->sub_xiter(); !jj.eof(); ++jj) {
		printf(" (row=%d : val=%g)", *jj, jj.val());
	}
	printf("\n");
}
@endcode

@see spsparse::VectorCooArray::dim_beginnings_xiter()
*/
template<class VectorCooArrayT>
class DimBeginningsXiter : public STLXiter<std::vector<size_t>::const_iterator>
{
public:
	SPSPARSE_LOCAL_TYPES(VectorCooArrayT);
	typedef std::vector<size_t>::const_iterator DimIterT;

protected:
	VectorCooArrayT const *arr;
	int index_dim;	// Dimension corresponding to our "rows"
	int val_dim;	// Dimension corresponding to our "cols"

public:
	/** @param _arr Matrix to iterate over (must live at least as long as this).
	@param _index_dim Dimension for the outer loop (eg: 0 for row major iteration).
	@param _val_dim Dimension for the inner loop (eg: 1 for row major iteration).
	@param dim_beginnings_begin Iterator in dim_beginnings array indicating the start of the outer loop (see spsparse::dim_beginnings()).
	@param dim_beginnings_end Iterator in dim_beginnings array indicating the end of the outer loop (see spsparse::dim_beginnings()). */
	DimBeginningsXiter(
		VectorCooArrayT const *_arr,
		int _index_dim, int _val_dim,
		DimIterT const &dim_beginnings_begin,
		DimIterT const &dim_beginnings_end);

	bool eof() { return ((ii+1) == end); }	// ii+1: Remember the sentinel at the end!

	index_type operator*()
		{ return arr->index(index_dim, *ii); }

	/** @brief Iterate along the current column (if we're scanning row
	major), or row (if we're scanning column major) of the matrix.

	@param _val_dim Dimension to report via operator*()  (1 for row major, 0 for column major).
	@see spsparse::DimIndexIter, spsparse::VectorCooArray::dim_begin()
	*/
	typedef ValSTLXiter<typename VectorCooArrayT::const_dim_iterator> sub_xiter_type;
	sub_xiter_type sub_xiter(int _val_dim = -1);

	// No val()
};

// -------------- Method Definitions
template<class VectorCooArrayT>
DimBeginningsXiter<VectorCooArrayT>::
	DimBeginningsXiter(
		VectorCooArrayT const *_arr,
		int _index_dim, int _val_dim,
		DimIterT const &dim_beginnings_begin,
		DimIterT const &dim_beginnings_end)
	: STLXiter<std::vector<size_t>::const_iterator>(dim_beginnings_begin, dim_beginnings_end),
		arr(_arr), index_dim(_index_dim), val_dim(_val_dim)
	{}

template<class VectorCooArrayT>
	typename DimBeginningsXiter<VectorCooArrayT>::sub_xiter_type DimBeginningsXiter<VectorCooArrayT>::sub_xiter(int _val_dim)
	{
		if (_val_dim < 0) _val_dim = val_dim;
		size_t a0 = *ii;
		size_t a1 = *(ii + 1);
		return sub_xiter_type(arr->dim_iter(_val_dim, a0), arr->dim_iter(_val_dim, a1));
	}


// -----------------------------------------------------
/** @brief Sorts an array and removes duplicates.
@param ret Accumulator for output.
@param A Input.
@param sort_order Order of dimensions to sort.  Use {0,1} for row major.
@param duplicate_policy What to do when duplicate entries are encountered (ADD (default), LEAVE_ALONE, REPLACE).
@param zero_nan If true, treat NaNs as zeros in the matrix (i.e. remove them).  This will prevent NaNs from propagating in computations.
*/
template<class VectorCooArrayT, class AccumulatorT>
void consolidate(AccumulatorT &ret,
	VectorCooArrayT const &A,
	std::array<int, VectorCooArrayT::rank> const &sort_order,
	DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
	bool zero_nan = false);	// Treat NaN like 0

template<class VectorCooArrayT, class AccumulatorT>
void consolidate(AccumulatorT &ret,
	VectorCooArrayT const &A,
	std::array<int, VectorCooArrayT::rank> const &sort_order,
	DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
	bool zero_nan = false)	// Treat NaN like 0
{
	const int RANK = VectorCooArrayT::rank;
	typedef typename VectorCooArrayT::index_type IndexT;
	typedef typename VectorCooArrayT::val_type ValT;

	// Nothing to do for zero-size matrices
	if (A.size() > 0) {

		// Get a sorted permutation
		std::vector<size_t> perm(sorted_permutation(A, sort_order));

		// Scan through to identify duplicates
		auto ii(perm.begin());

		// Skip over initial 0 and NaN
		for (;; ++ii) {
			if (ii == perm.end()) goto finished;		// Nothing more to do
			if (!isnone(A.val(*ii), zero_nan)) break;
		}

		// New temporary entry
		std::array<IndexT,RANK> accum_idx(A.index(*ii));
		ValT accum_val = A.val(*ii);
		++ii;

		for (; ; ++ii) {
			// Skip over initial 0 and NaN
			for (;; ++ii) {
				if (ii == perm.end()) {
					// Write out the last thing we had in our accumulator
					ret.add(accum_idx, accum_val);

					goto finished;		// Nothing more to do
				}
				if (!isnone(A.val(*ii))) break;
			}

			// Test if A.index(*ii) == accum_idx
			auto new_idx = A.index(*ii);
			for (int k=0; k<RANK; ++k) {
				if (new_idx[k] != accum_idx[k]) {
					// They don't match.  Write out our accumulator, and reset to this one
					ret.add(accum_idx, accum_val);
					accum_idx = new_idx;
					accum_val = A.val(*ii);
					goto continue_outer;
				}
			}

			// They do match!  Add it into the accumulator...
			if (duplicate_policy == DuplicatePolicy::ADD)
				accum_val += A.val(*ii);
			else if (duplicate_policy == DuplicatePolicy::REPLACE)
				accum_val = A.val(*ii);

		continue_outer: ;
		}
	}
finished:


	ret.set_sorted(sort_order);
}
// -----------------------------------------------------
/** @brief Consolidates an array, if it is not already consolidated.
Does not overwrite the original.
@see spsparse::consolidate() */
template<class ArrayT>
class Consolidate
{
	static const int rank = ArrayT::rank;
	std::unique_ptr<ArrayT> A2;		// If needed.
	ArrayT const *Ap;						// The final answer
public:
	/** 
	@param A Matrix to consolidate. (must live at least as long as this).
	@param sort_order Order of dimensions to sort.  Use {0,1} for row major.
	@param duplicate_policy What to do when duplicate entries are encountered (ADD (default), LEAVE_ALONE, REPLACE).
	@param zero_nan If true, treat NaNs as zeros in the matrix (i.e. remove them).  This will prevent NaNs from propagating in computations.
	@see spsparse::consolidate() */
	Consolidate(ArrayT const *A,
		std::array<int, ArrayT::rank> const &sort_order,
		DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
		bool zero_nan = false);

	/** @brief Produces the consolidated array.

	This might be the original array called in the constructor, if
	that was properly consolidated; or it might be a new array created
	in this class. */
	ArrayT const &operator()() {
		return *Ap;
	}
};

// -------------------- Method Definitions
template<class ArrayT>
Consolidate<ArrayT>::
	Consolidate(ArrayT const *A,
		std::array<int, ArrayT::rank> const &sort_order,
		DuplicatePolicy duplicate_policy,
		bool zero_nan)
	{
		if (A->sort_order == sort_order) {
			// A is already consolidated, use it...
			Ap = A;
		} else {
			// Must consolidate A...
			A2 = A->new_blank();
			Ap = A2.get();
			consolidate(*A2, *A, sort_order, duplicate_policy, zero_nan);
		}
	}


// -------------------------------------------------------------
// --------------------------------------------------------
/** @brief Internal class used by spsparse::sorted_permutation(). */
template<class VectorCooArrayT>
struct CmpIndex {
	const int RANK = VectorCooArrayT::rank;

	VectorCooArrayT const *arr;
	std::array<int, VectorCooArrayT::rank> sort_order;

	CmpIndex(VectorCooArrayT const *_arr,
		std::array<int, VectorCooArrayT::rank> const &_sort_order)
	: arr(_arr), sort_order(_sort_order) {}

	bool operator()(int i, int j)
	{
		for (int k=0; k<RANK-1; ++k) {
			int dim = sort_order[k];
			if (arr->index(dim,i) < arr->index(dim,j)) return true;
			if (arr->index(dim,i) > arr->index(dim,j)) return false;
		}
		int dim = sort_order[RANK-1];
		return arr->index(dim,i) < arr->index(dim,j);
	}
};
// --------------------------------------------------------------------
/** @brief Generates a permutation that, if applied, would result in the array being sorted.

@param A Input array.
@param sort_order Order of dimensions to sort.  Use {0,1} for row major.
@return The permutation.

@note Sorting is done in-place.  Elements added first will remain
      first, allowing for proper behavior of duplicate_policy in
      spsparse::consolidate(). */
template<class VectorCooArrayT>
std::vector<size_t> sorted_permutation(VectorCooArrayT const &A,
	std::array<int, VectorCooArrayT::rank> const &sort_order);

template<class VectorCooArrayT>
std::vector<size_t> sorted_permutation(VectorCooArrayT const &A,
	std::array<int, VectorCooArrayT::rank> const &sort_order)
{
	// Decide on how we'll sort
	CmpIndex<VectorCooArrayT> cmp(&A, sort_order);

	// Generate a permuatation
	int n = A.size();
	std::vector<size_t> perm; perm.reserve(n);
	for (int i=0; i<n; ++i) perm.push_back(i);

	// Sort the permutation
	std::stable_sort(perm.begin(), perm.end(), cmp);

	return perm;
}
// ----------------------------------------------------------
template<class TypeT, int RANK, class CooArrayT>
void to_sparse(CooArrayT &ret,
	blitz::Array<TypeT, RANK> const &arr);

template<class TypeT, int RANK, class CooArrayT>
void to_sparse(CooArrayT &ret,
	blitz::Array<TypeT, RANK> const &arr)
{
	for (auto ii=arr.begin(); ii != arr.end(); ++ii) {
		if (*ii != 0) ret.add_blitz(ii.position(), *ii);
	}
}
// ----------------------------------------------------------

/** @} */
}	// Namespace

#endif // Guard
