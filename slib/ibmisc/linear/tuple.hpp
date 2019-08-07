#ifndef IBMISC_LINEAR_TUPLE_HPP
#define IBMISC_LINEAR_TUPLE_HPP

#include <memory>
#include <Eigen/SparseCore>
#include <spsparse/SparseSet.hpp>
#include <ibmisc/linear/linear.hpp>
#include <spsparse/eigen.hpp>

using namespace spsparse;

namespace ibmisc {
namespace linear {


/** Return value of a sparse matrix */
struct Weighted_Tuple : public Weighted {

    template<int RANK>
        using TupleListLT = spsparse::TupleList<long,double,RANK>;

    TupleListLT<1> wM;
    TupleListLT<2> M;
    TupleListLT<1> Mw;


    /** Construct with shared dimensions from elsewhere */
    Weighted_Tuple(bool conservative=true)
        : Weighted(LinearType::TUPLE, conservative) {}

    void set_shape(std::array<long,2> _shape);

    /** Sparse shape of the matrix */
    std::array<long,2> shape() const
        { return M.shape(); }

    /** Takes inner product between a weight vector and the input(s) As.
    @param dim 0=B (output) weight, 1=A (intput) weight.
    @param As As[nvec,nA] Vector(s) to compute inner product.
    @param out out[nvec] Place output here.
    @param zero_out If set, zero indices in out that we will change, before applying.
        Indices for which weight==0 will not be touched. */
    void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true) const;

    /** Compute M * As.  (Sparse indexing)
    Does not touch the nullspace of M.
    NOTE: For in-place multiplication, see apply_M_inplace(). */
    void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) const;

    /** Computes As = M * As, result put in place. */
    virtual void apply_M_inplace(
        blitz::Array<double,2> &As,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) const;

    long nnz() const { return M.size(); }

protected:
    virtual void _to_coo(
        blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
        blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
        blitz::Array<double,1> &values) const;      // Must bepre-allocated(nnz)


    virtual void _get_weights(
        int idim,    // 0=wM, 1=Mw
        blitz::Array<double,1> &w) const;

    /** I/O */
    virtual void ncio(NcIO &ncio, std::string const &vname); 
public:
};


}}    // namespace
#endif
