#ifndef IBMISC_MATH_HPP
#define IBMISC_MATH_HPP

#include <cmath>
#include <limits>

namespace ibmisc {

/** Rounds a specific number of binary digits off the mantissa of a
floating-point number.  NOTES:

1. This does not twiddle IEEE bits directly.  It is compliant with the
   C++11 standard.

2. This function ROUNDS, it does not do a floor or ceiling.

@param x The number to round

@param ntruncate The number of low-order bits in the mantiss that will
be set to zero.  For example... an IEEE double has 53 bits in the
mantissa.  Setting ntruncat=3 will result in a number with 50 bits of
precision.
*/
inline double round_mantissa(double const x, int const ntruncate=3)
{
    int const digits = std::numeric_limits<double>::digits;
    int exp;
    double const mantissa = frexp(x, &exp);
    return ldexp(std::round(ldexp(x, digits-ntruncate-exp)),  +exp+ntruncate-digits);
}

}    // namespace
#endif
