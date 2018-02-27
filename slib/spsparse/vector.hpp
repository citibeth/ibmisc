#ifndef SPSPARSE_VECTOR_HPP
#define SPSPARSE_VECTOR_HPP

#include <vector>

namespace spsparse {
namespace vaccum {

template<class ValueT>
class Vector {
public:
    typedef ValueT val_type;
private:
    std::vector<ValueT> &vec;
public:
    Vector(std::vector<ValueT> &_vec) : vec(_vec) {}
    void add(ValueT const &val)
        { vec.push_back(val); }
};

template<class ValueT>
inline Vector<ValueT> vector(std::vector<ValueT> &vals)
    { return Vector<ValueT>(vals); }

}}    // namespace

#endif
