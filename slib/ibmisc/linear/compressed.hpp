#ifndef IBMISC_LINEAR_COMPRESSED_HPP
#define IBMISC_LINEAR_COMPRESSED_HPP

#include <ibmisc/linear/linear.hpp>
#include <ibmisc/zarray.hpp>

namespace ibmisc {
namespace linear {

class Weighted_Eigen;

// ==================================================================
class Weighted_Compressed : public Weighted
{
    friend Weighted_Compressed compress(Weighted_Eigen &eigen);

    std::array<ZArray<int,double,1>, 2> weights;    // {wM, Mw}
    ZArray<int,double,2> M;

public:
    Weighted_Compressed() : Weighted(LinearType::COMPRESSED) {}

    // ================= Implements Weighted
    /** Sparse shape of the matrix */
    std::array<long,2> shape() { return M.shape(); }

    /** Computes out = As * weights[dim] */
    void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true);

    /** Computes out = M * As */
    void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true);

    void ncio(NcIO &ncio, std::string const &vname);

};

extern Weighted_Compressed compress(Weighted_Eigen &eigen);




}};    // namespace
#endif    // guad