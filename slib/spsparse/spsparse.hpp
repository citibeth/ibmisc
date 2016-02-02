#ifndef SPSPARSE_SPSPARSE_H
#define SPSPARSE_SPSPARSE_H

#include <exception>
#include <iostream>

#include <ibmisc/blitz.hpp>

namespace spsparse {

/** @defgroup spsparse spsparse.hpp
@brief Basic stuff common to all of SpSparse library.

@{
*/


/** @brief What to do in algorithms when duplicate entries are encountered.

- ADD (default): Sum them together.
- LEAVE_ALONE: Use the first value encountered.
- REPLACE: Use the last value encountered.

@see spsparse::consolidate() */
enum class DuplicatePolicy {
	LEAVE_ALONE, ADD, REPLACE};


/** @brief Excpetion thrown by the default SpSparse error handler. */
class Exception : public std::exception
{
public:
	virtual ~Exception()
		{}

	virtual const char* what() const noexcept
		{ return "spsparse::Exception()"; }
};


// This is not possible in C++11
//extern std::function<void(int, const char *, ...)> error;

// Use this instead.
// http://www.thecodingforums.com/threads/function-pointers-to-printf.317925/
/** @brief Printf-like signature of error handle functions to be used by SpSparse. */
typedef void (*error_ptr) (int retcode, char const *str, ...);

/** @brief Error handler used by SpSparse.  May be changed by user's
main program, to fit into some larger error handling system (eg:
Everytrace).

https://github.com/citibob/everytrace */
extern error_ptr spsparse_error;

/** @brief Promote relevant template parameters.

Used to propagate relevant template parameters throughout the
different classes that need them.  Provides: rank, index_type, val_type, indices_type and blitz_type.

@note val_type is DIFFERENT from the standard STL value_type.  Standard STL value_type is the same as our index_type.

Code Example
@code
template<class VectorCooArrayT>
class MyClass {
public:
	SPSPARSE_LOCAL_TYPES(VectorCooArrayT);

};
@endcode
*/
#define SPSPARSE_LOCAL_TYPES(ArrayOrIterT) \
	static const int rank = ArrayOrIterT::rank; \
	typedef typename ArrayOrIterT::index_type index_type; \
	typedef typename ArrayOrIterT::val_type val_type; \
	typedef std::array<index_type, rank> indices_type; \
	typedef blitz::Array<val_type, rank> blitz_type

// -------------------------------------------------------------
// Values for sort_order formal parameter below
extern const std::array<int,2> ROW_MAJOR;
extern const std::array<int,2> COL_MAJOR;

// The expression "std::isnan(n) || (n == 0)" for different data types.
// Use template specilization here...
/** @brief Internal helper function.

Used to tell if a value is "none" (i.e. 0 or NaN) and should therefore
be eliminated under zero_nan.  This will need to be specialized for
value types that are not a single IEEE floating point (eg:
std::complex<double>).

@see spsparse::consolidate() */
template<class NumT>
inline bool isnone(NumT const n, bool const zero_nan=false)
{
	if (zero_nan) {
		return std::isnan(n) || (n == 0);
	} else {
		return (n == 0);
	}
}

// ----------------------------------------------------------
} // Namespace

// -------------------------------------------------------------
/** Hack to write std::array to ostream. */
template<class T>
std::ostream &stream(std::ostream &os, T const * const a, int RANK);

template<class T>
std::ostream &stream(std::ostream &os, T const * const a, int RANK)
{
	if (RANK == 0) {
		os << "{}";
	} else {
		os << "{";
		for (int k=0; ; ) {
			os << a[k];
			++k;
			if (k == RANK) break;
			os << ", ";
		}
		os << "}";
	}
	return os;
}

/** @} */

#endif // Guard
