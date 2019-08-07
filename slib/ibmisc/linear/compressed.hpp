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

public:
    std::array<ZArray<int,double,1>, 2> weights;    // {wM, Mw}
    ZArray<int,double,2> M;

public:
    Weighted_Compressed() : Weighted(LinearType::COMPRESSED) {}

    void set_shape(std::array<long,2> _shape);

    // ================= Implements Weighted
    /** Sparse shape of the matrix */
    std::array<long,2> shape() const;

    /** Computes out = As * weights[dim] */
    void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true) const;

    /** Computes out = M * As
    NOTE: As and out cannot be the same! */
    void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) const;

    void ncio(NcIO &ncio, std::string const &vname);

    long nnz() const;

protected:
    void _to_coo(
        blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
        blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
        blitz::Array<double,1> &values) const;      // Must bepre-allocated(nnz)

    void _get_weights(
        int idim,    // 0=wM, 1=Mw
        blitz::Array<double,1> &w) const;
};

extern Weighted_Compressed compress(Weighted_Eigen &eigen);







}};    // namespace
#endif    // guad
