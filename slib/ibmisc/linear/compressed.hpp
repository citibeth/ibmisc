#ifndef IBMISC_LINEAR_COMPRESSED_HPP
#define IBMISC_LINEAR_COMPRESSED_HPP

#include <ibmisc/linear/linear.hpp>
#include <ibmisc/zarray.hpp>

namespace ibmisc {
namespace linear {

// ==================================================================
class Weighted_Compressed : public Weighted
{
    std::array<ZArray<int,double,1>, 2> weights;    // {wM, Mw}
    ZArray<int,double,2> M;

    Weighted_Compressed() : Weighted(LinearType::COMPRESSED) {}

    // ================= Implements Weighted
protected:
    /** Computes out = As * weights[dim] */
    void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true);

public:

    /** Computes out = M * As */
    void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        FillType fill_type=FillType::nan_all,
        bool force_conservation=true);

    void ncio(NcIO &ncio, std::string const &vname);

};


}};    // namespace
#endif    // guad
