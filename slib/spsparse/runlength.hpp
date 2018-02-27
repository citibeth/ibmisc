#ifndef IBMISC_RUNLENGTH_HPP
#define IBMISC_RUNLENGTH_HPP

#include <cmath>
#include <blitz/array.h>
#include <boost/enum.hpp>

namespace spsparse {

BOOST_ENUM_VALUES(RLAlgo
    (PLAIN) (0)
    (DIFFS) (1)
)

// -----------------------------------------------------------------
/** Use std::equal_to<TypeT> instead, if you don't like this. */
template<class TypeT>
struct DefaultRLEqual : public std::equal_to<TypeT>
{}

template<>
void DefaultRLEqual<double>::operator()(double const &a, double const &b) const
{
    switch(
        (std::isnan(a) ? 2 : 0) +
        (std::isnan(b) ? 1 : 0))
    {
        case 3:
            return true;
        case 1:
        case 2:
            return false;
        case 0:
            return (a == b);
    }
}
// -----------------------------------------------------------------
// -------------------------------------------------------------
namespace vaccum {
/** Runlength encodes to a Scalar accumulator (uses just add(x)) */
template <class CountsVAccumT, class ValuesVAccumT, class EqualT>
class RLEncode
{
    typedef typename CountsVAccumT::val_type count_type;
    typedef typename ValuesVAccumT::val_type value_type;
public:
    CountsVAccumT counts_vaccum;
    ValuesVAccumT values_vaccum;
    RLAlgo algo;
    EqualT eq;

    value_type last_raw;    // Used only for difference encoding
    value_type run_val;
    bool first = true;
    count_type count;

    RLEncode(
        CountsVAccumT &&_counts_vaccum,
        ValuesVAccumT &&_values_vaccum,
        RLAlgo _algo,
        EqualT const &&_eq)
    : counts_vaccum(std::move(_counts_vaccum)),
      values_vaccum(std::move(_values_vaccum)),
      algo(_algo), eq(std::move(_eq)) {}

    void add(value_type raw)
    {
        if (first) {
            run_val = raw;
            last_raw = raw;
            count = 1;
            first = false;
        } else {

            // Do difference encoding if requested
            value_type val;
            switch(algo.index()) {
                case RLAlgo::PLAIN :
                    val = raw - last_raw;
                    last_raw = raw;
                break;
                case RLAlgo::DIFFS :
                    val = raw;
                break;
            }

            // Runlength Encoding
            if (eq(val, run_val)) {
                ++count;
            } else {
                counts_vaccum.add(count);
                values_vaccum.add(run_val);
                run_val = val;
                count = 1;
            }
        }
    }

    ~RLEncode()
    {
        if (first) return;
        counts_vaccum.add(count);
        values_vaccum.add(run_val);
    }

};

template <class CountsVAccumT, class ValuesVAccumT, class EqualT=DefaultRLEqual<typename ValuesVAccumT::val_type>>
inline RLEncode<CountsVAccumT, ValuesVAccumT, EqualT>
rl_encode(
    CountsVAccumT &&counts_vaccum,
    ValuesVAccumT &&values_vaccum,
    RLAlgo algo = RLAlgo::PLAIN,
    EqualT const &&eq=std::equal_to<typename ValuesVAccumT::val_type>())
{
    return RLEncode<CountsVAccumT, ValuesVAccumT, EqualT>(
        std::move(counts_vaccum),
        std::move(values_vaccum),
        algo, std::move(eq));
}
}    // namespace vaccum
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------



// =============================================================
/** Usage:
       for (RLDecode rl(...); ++rl; ) print(*rl);
*/
template<class CountsIterT, class ValuesIterT>
class RLDecode
{
    typedef typename ValuesIterT::value_type value_type;
    RLAlgo const algo;

    value_type cur_raw = 0;
    value_type cur_val;
    typename CountsIterT::value_type cur_count;

    CountsIterT counts_iter, counts_end;
    ValuesIterT values_iter, values_end;
public:
    RLDecode(
        CountsIterT &&_counts_begin, CountsIterT &&_counts_end,
        ValuesIterT &&_values_begin, ValuesIterT &&_values_end,
        RLAlgo _algo)
    : algo(_algo), cur_raw(0), cur_val(0), cur_count(0),
        counts_iter(std::move(_counts_begin)),
        counts_end(std::move(_counts_end)),
        values_iter(std::move(_values_begin)),
        values_end(std::move(_values_end))
    {}

    bool operator++()
    {
        for (;;) {
            if (cur_count > 0) {
                --cur_count;
                switch(algo) {
                    case RLAlgo::DIFFS :
                        cur_raw += cur_val;
                    break;
                }
                return true;    // Have a value; get it with operator*()
            } else {
                if (counts_iter == counts_end) return false;    // No more values
                cur_count = *counts_iter;
                cur_val = *values_iter;
                ++counts_iter;
                ++values_iter;
            }
        }
    }

    value_type const &operator*() const
    {
        switch(algo) {
            case RLAlgo::PLAIN :
                return cur_plain;
            case RLAlgo::DIFFS:
                return cur_raw;
        }
    }
};


template<class CountsIterT, class ValuesIterT>
RLDecode<CountsIterT, ValuesIterT> rl_decode(
    CountsIterT &&counts_begin, CountsIterT &&counts_end,
    ValuesIterT &&values_begin, ValuesIterT &&values_end,
    RLAlgo algo = RLAlgo::PLAIN)
{
    return RLDecode<CountsIterT, ValuesIterT>(
        std::move(counts_begin), std::move(counts_end),
        std::move(values_begin), std::move(values_end),
        algo);
}

// ===========================================================================





// ---------------------------------------------------------------------
} // namespace spsparse
#endif
