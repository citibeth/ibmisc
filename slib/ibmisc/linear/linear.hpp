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


class Weighted {
public:
    // TmpAlloc tmp;    // Sometimes, hold the things we're wrapping.
    LinearType const type;

    /** True if this regridding matrix is conservative.  Matrices could be
    non-conservative, for example, in the face of smoothing on I.  Or when
    regridding between the IceBin and ModelE ice sheets. */
    bool conservative;    // Is this matrix conservative?

protected:
    Weighted(LinearType _type) : type(_type), conservative(true) {}
    Weighted(LinearType _type, bool _conservative) : type(_type), conservative(_conservative) {}

    virtual ~Weighted() {}

    virtual void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true) = 0;

public:
    /** Sparse shape of the matrix */
    virtual std::array<long,2> shape() = 0;

    virtual void ncio(NcIO &ncio, std::string const &vname); 

    /** Compute M * As.
    Does not touch the nullspace of M. */
    virtual void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) = 0;

    /** Compute wM * As */
    void apply_wM(
        blitz::Array<double,2> const &As,
        blitz::Array<double,1> &out,
        AccumType accum_type=AccumType::REPLACE)
    { apply_weight(0, As, out, accum_type, false); }

    /** Compute Mw * As */
    void apply_Mw(
        blitz::Array<double,2> const &As,
        blitz::Array<double,1> &out,
        AccumType accum_type=AccumType::REPLACE)
    { apply_weight(1, As, out, accum_type, false); }

    /** Compute 1. / (wM * As) */
    void apply_sM(
        blitz::Array<double,2> const &As,
        blitz::Array<double,1> &out,
        AccumType accum_type=AccumType::REPLACE)
    { apply_weight(0, As, out, accum_type, true); }

    /** Compute 1. / (Mw * As) */
    void apply_Ms(
        blitz::Array<double,2> const &As,
        blitz::Array<double,1> &out,
        AccumType accum_type=AccumType::REPLACE)
    { apply_weight(1, As, out, accum_type, true); }

};

extern std::unique_ptr<Weighted> new_weighted(LinearType type);

std::unique_ptr<Weighted> nc_read_weighted(netCDF::NcGroup *nc, std::string const &vname);


}};    // namespace
#endif    // guad
