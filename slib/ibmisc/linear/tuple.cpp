#include <ibmisc/linear/tuple.hpp>

namespace ibmisc {
namespace linear {

void Weighted_Tuple::apply_weight(
    int dim,    // 0=B, 1=A
    blitz::Array<double,2> const &As,    // As(nvec, ndim)
    blitz::Array<double,1> &out,
    bool zero_out) const
{
    (*ibmisc_error)(-1, "linear::Weighted_Tuple::apply_weight() not yet implemented");
}

void Weighted_Tuple::apply_M(
    blitz::Array<double,2> const &As,
    blitz::Array<double,2> &out,
    AccumType accum_type,
    bool force_conservation) const
{
    (*ibmisc_error)(-1, "linear::Weighted_Tuple::apply_M() not yet implemented");
}

void Weighted_Tuple::apply_M_inplace(
    blitz::Array<double,2> &As,
    AccumType accum_type,
    bool force_conservation) const
{
    (*ibmisc_error)(-1, "linear::Weighted_Tuple::apply_M_inplace() not yet implemented");
}

void Weighted_Tuple::_to_coo(
    blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
    blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
    blitz::Array<double,1> &values) const      // Must bepre-allocated(nnz)
{
    (*ibmisc_error)(-1, "linear::Weighted_Tuple::to_coo() not yet implemented");
}

void Weighted_Tuple::_get_weights(
    int idim,    // 0=wM, 1=Mw
    blitz::Array<double,1> &w) const
{
    (*ibmisc_error)(-1, "linear::Weighted_Tuple::get_weights() not yet implemented");
}

/** I/O */
void Weighted_Tuple::ncio(NcIO &ncio, std::string const &vname)
{
}

}}    // namespace
