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

} // namespace
#endif
