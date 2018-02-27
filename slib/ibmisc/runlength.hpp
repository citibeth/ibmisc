#ifndef IBMISC_RUNLENGTH_HPP
#define IBMISC_RUNLENGTH_HPP

#include <cmath>
#include <blitz/array.h>

namespace ibmisc {

template<class TypeT>
struct RLEncode {
    std::vector<int> ends;
    std::vector<TypeT> values;
};

struct EqualUsingNaN {
    bool operator()(double a, double b) const {
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
};

template<class TypeT,class EqualT=std::equal_to<TypeT>>
RLEncode<TypeT> rlencode(
    blitz::Array<TypeT,1> const &vec,
    bool diff_encode = false,
    EqualT const &eq=std::equal_to<TypeT>());

template<class TypeT,class EqualT>
RLEncode<TypeT> rlencode(
    blitz::Array<TypeT,1> const &vec,
    bool diff_encode,
    EqualT const &eq)
{
    RLEncode<TypeT> ret;

    if (vec.extent(0) == 0) return ret;

    TypeT run_val = vec(0);
    int i=1;
    for (; i<vec.extent(0); ++i) {
        TypeT const vec_i = (diff_encode ? vec(i)-vec(i-1) : vec(i));
//std::cout << " " << vec_i;
        if (!(eq(vec_i, run_val))) {
            ret.ends.push_back(i);
            ret.values.push_back(run_val);
            run_val = vec_i;
        }
    }
//std::cout << std::endl;
    ret.ends.push_back(i);
    ret.values.push_back(run_val);

    return ret;
}

template<class TypeT>
blitz::Array<TypeT,1> rldecode(
    blitz::Array<int,1> const &ends,
    blitz::Array<TypeT,1> const &values,
    bool diff_encode = false)
{
    if (ends.extent(0) == 0) return blitz::Array<TypeT,1>(0);

    blitz::Array<TypeT,1> ret(ends(ends.size()-1));
    int j=0;
    TypeT last_ret_j = 0;
    for (int i=0; i<ends.size(); ++i) {
        for (; j<ends(i); ++j) {
            if (diff_encode) {
                last_ret_j += values(i);
                ret(j) = last_ret_j;
            } else {
                ret(j) = values(i);
            }
        }
    }

    return ret;

}

// ----------------------------------------------------------
// -------------------------------------------------------------

/** Runlength encodes to a Scalar accumulator (uses just add(x)) */
template <class CountsAccumT, class ValuesAccumT, class EqualT>
class RLEncoder
{
    typedef typename CountsAccumT::val_type count_type;
    typedef typename ValuesAccumT::val_type value_type;
public:
    CountsAccumT counts_accum;
    ValuesAccumT values_accum;
    bool diff_encode;
    EqualT eq;

    value_type last_raw;    // Used only for difference encoding
    value_type run_val;
    bool first = true;
    count_type count;

    RLEncoder(
        CountsAccumT &&_counts_accum,
        ValuesAccumT &&_values_accum,
        bool _diff_encode,
        EqualT const &&_eq)
    : counts_accum(std::move(_counts_accum)),
      values_accum(std::move(_values_accum)),
      diff_encode(_diff_encode), eq(std::move(_eq)) {}

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
            if (diff_encode) {
                val = raw - last_raw;
                last_raw = raw;
            } else {
                val = raw;
            }

            // Runlength Encoding
            if (eq(val, run_val)) {
                ++count;
            } else {
                counts_accum.add(count);
                values_accum.add(run_val);
                run_val = val;
                count = 1;
            }
        }
    }

    ~RLEncoder()
    {
        if (first) return;
        counts_accum.add(count);
        values_accum.add(run_val);
    }

};

template <class CountsAccumT, class ValuesAccumT, class EqualT=std::equal_to<typename ValuesAccumT::val_type>>
RLEncoder<CountsAccumT, ValuesAccumT, EqualT>
inline rl_encoder(
    CountsAccumT &&counts_accum,
    ValuesAccumT &&values_accum,
    bool diff_encode = false,
    EqualT const &&eq=std::equal_to<typename ValuesAccumT::val_type>())
{
    return RLEncoder<CountsAccumT, ValuesAccumT, EqualT>(
        std::move(counts_accum),
        std::move(values_accum),
        diff_encode, std::move(eq));
}



// =============================================================
/** Usage:
       for (RLDecoder rl(...); ++rl; ) print(*rl);
*/
template<class CountsIterT, class ValuesIterT>
class RLDecoder
{
    typedef typename ValuesIterT::value_type value_type;
    bool const diff_encode;

    value_type cur_raw = 0;
    value_type cur_val;
    typename CountsIterT::value_type cur_count;

    CountsIterT counts_iter, counts_end;
    ValuesIterT values_iter, values_end;
public:
    RLDecoder(
        CountsIterT &&_counts_begin, CountsIterT &&_counts_end,
        ValuesIterT &&_values_begin, ValuesIterT &&_values_end,
        bool _diff_encode)
    : diff_encode(_diff_encode), cur_raw(0), cur_val(0), cur_count(0),
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
                if (diff_encode) cur_raw += cur_val;
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
        return (diff_encode ? cur_raw : cur_val);
    }
};


template<class ValuesIterT, class CountsIterT>

RLDecoder<CountsIterT, ValuesIterT> rl_decoder(
    CountsIterT &&counts_begin, CountsIterT &&counts_end,
    ValuesIterT &&values_begin, ValuesIterT &&values_end,
    bool diff_encode = false)
{
    return RLDecoder<CountsIterT, ValuesIterT>(
        std::move(counts_begin), std::move(counts_end),
        std::move(values_begin), std::move(values_end),
        diff_encode);
}


} // namespace
#endif
