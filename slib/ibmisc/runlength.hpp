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
template<class ValueT>
class VectorAccum {
    std::vector<ValueT> &vec;
    VectorAccum(std::vector<ValueT> &_vec) : vec(_vec) {}
    void add(std::vector<ValueT> const &val)
        { vec.push_back(val); }
};

template<class ValueT>
inline VectorAccum<ValueT> vector_accum(std::vector<ValueT> const &val)
    { return VectorAccum<ValueT>(val); }
// -------------------------------------------------------------

/** Runlength encodes to a Scalar accumulator (uses just add(x)) */
template <class ValuesAccumT, class CountsAccumT, class EqualT=std::equal_to<typename ValuesAccumT::value_type>>
class RLEncoder
{
    typedef typename CountsAccumT::value_type count_type;
    typedef typename ValuesAccumT::value_type value_type;
public:
    ValuesAccumT values_accum;
    CountsAccumT counts_accum;
    bool diff_encode;
    EqualT eq;

    value_type last_raw;    // Used only for difference encoding
    value_type run_val;
    bool first = true;
    count_type count;

    RLEncoder(
        ValuesAccumT &&_values_accum,
        CountsAccumT &&_counts_accum,
        bool _diff_encode = false,
        EqualT const &&_eq=std::equal_to<typename ValuesAccumT::value_type>())
    : values_accum(std::move(_values_accum)),
        counts_accum(std::move(_counts_accum)),
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
            if (!eq(val, run_val)) {
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

template <class ValuesAccumT, class CountsAccumT, class EqualT=std::equal_to<typename ValuesAccumT::value_type>>
RLEncoder<ValuesAccumT, CountsAccumT, EqualT>
inline rl_encoder(
    ValuesAccumT &&values_accum,
    CountsAccumT &&counts_accum,
    bool diff_encode = false,
    EqualT const &&eq=std::equal_to<typename ValuesAccumT::value_type>())
{
    return RLEncoder<ValuesAccumT, CountsAccumT, EqualT>(
        std::move(values_accum), std::move(counts_accum),
        diff_encode, std::move(eq));
}



// =============================================================
/** Usage:
       for (RLGenerator rl(...); ++rl; ) print(*rl);
*/
template<class ValuesIterT, class CountsIterT>
class RLGenerator
{
    typedef typename ValuesIterT::value_type value_type;
    bool const diff_encode;

    value_type cur_raw = 0;
    value_type cur_val;
    typename CountsIterT::value_type cur_count;

    ValuesIterT values_iter, values_end;
    CountsIterT counts_iter, counts_end;
public:
    RLGenerator(ValuesIterT &&_values_begin, ValuesIterT &&_values_end,
        CountsIterT &&_counts_begin, CountsIterT &&_counts_end,
        bool _diff_encode = false)
    : diff_encode(_diff_encode), cur_raw(0), cur_val(0), cur_count(0),
    values_iter(std::move(_values_begin)), counts_iter(std::move(_counts_begin)),
    values_end(std::move(_values_end)), counts_end(std::move(_counts_end))
    {
        ++counts_iter;
        ++values_iter;
    }

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





} // namespace
#endif
