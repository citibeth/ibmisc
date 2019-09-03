#ifndef IBMISC_LINEAR_LINEAR_HPP
#define IBMISC_LINEAR_LINEAR_HPP

#include <boost/enum.hpp>
#include <blitz/array.h>
#include <ibmisc/netcdf.hpp>

namespace ibmisc {
namespace linear {

BOOST_ENUM_VALUES(LinearType, int,
    (EIGEN) (0)
    (COMPRESSED) (1)
    (TUPLE) (2)
)

// What do do with output values in the active space
BOOST_ENUM_VALUES(AccumType, int,
    (REPLACE) (0)
    (ACCUMULATE) (1)        // Accumulate; error if any NaNs creep in
    (REPLACE_OR_ACCUMULATE) (2)    // Replace if NaN, else accumulate
)



#if 0
/** Abstract sparse vector, just enough to multiply by. */
class Vector {
    virtual ~Vector() {}

    /** Compute inner product with v
    @param v_s Other vector to compute inner product with.
        Indexing!
    @param invert If true, then compute 1/inner product. */
    virtual void apply(
        blitz::Array<double,2> const &As,
        blitz::Array<double,1> &out,
        FillType accum_type,
        bool invert = false) = 0;

};

/** Abstract "compressed" matrix, just enough to multiply by. */
class Matrix
{
    virtual ~Matrix() {}

    /** Computes M*v_s, stores result in out_s. */
    virtual void apply(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        FillType accum_type) = 0;

};
#endif

/** Encapsulates a (possibly not conservative) regridding matrix with
weight vectors for input and output vector space.  Implementations use
either with Eigen matrices on a subspace, or with Zlib-compressed
matrices. */
class Weighted {
public:
    // TmpAlloc tmp;    // Sometimes, hold the things we're wrapping.
    LinearType type;

    /** True if this regridding matrix is conservative.  Matrices could be
    non-conservative, for example, in the face of smoothing on I.  Or when
    regridding between the IceBin and ModelE ice sheets. */
    bool conservative;    // Is this matrix conservative?

    /** True if this matrix is scaled */ 
    bool scaled;

protected:
    Weighted(LinearType _type, bool _conservative=true, bool _scaled=false)
        : type(_type), conservative(_conservative), scaled(_scaled) {}

public:
    virtual ~Weighted() {}

    /** Sparse shape of the matrix */
    virtual std::array<long,2> shape() const = 0;

    /** Takes inner product between a weight vector and the input(s) As.
    @param dim 0=B (output) weight, 1=A (intput) weight.
    @param As As[nvec,nA] Vector(s) to compute inner product.
    @param out out[nvec] Place output here.
    @param zero_out If set, zero indices in out that we will change, before applying.
        Indices for which weight==0 will not be touched. */
    virtual void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true) const = 0;

    /** Compute M * As.  (Sparse indexing)
    Does not touch the nullspace of M.
    NOTE: For in-place multiplication, see apply_M_inplace(). */
    virtual void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) const = 0;

    /** Computes As = M * As, result put in place. */
    virtual void apply_M_inplace(
        blitz::Array<double,2> const &As,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) const
    {
        (*ibmisc_error)(-1,
            "apply_M_inplace() not implemented for this linear::Weighted type!");
    }


    /** @return {weights[0].nnz, M.nnz, weights[1].nnz} */
    virtual long nnz() const = 0;

protected:
    virtual void _to_coo(
        blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
        blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
        blitz::Array<double,1> &values) const = 0;      // Must bepre-allocated(nnz)


    virtual void _get_weights(
        int idim,    // 0=wM, 1=Mw
        blitz::Array<double,1> &w) const = 0;

public:

//    virtual void map(std::function<void(std::array<int,2>, double)> const &fn) const = 0;

    /** @return The sparse matrix in uncompressed, easily convertible form. */
    void to_coo(
        blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
        blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
        blitz::Array<double,1> &values) const;

    void get_weights(
        int idim,    // 0=wM, 1=Mw
        blitz::Array<double,1> &w) const;

    /** I/O */
    virtual void ncio(NcIO &ncio, std::string const &vname); 
};

extern std::unique_ptr<Weighted> new_weighted(LinearType type);

/** Read a Weighted matrix, either eigen OR compressed. */
std::unique_ptr<Weighted> nc_read_weighted(netCDF::NcGroup *nc, std::string const &vname);


}};    // namespace
#endif    // guad
