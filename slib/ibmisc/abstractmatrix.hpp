#ifndef IBMISC_SPARSEMATRIX_HPP
#define IBMISC_SPARSEMATRIX_HPP

namespace ibmisc {

#if 0
/** Abstract sparse vector, just enough to multiply by. */
template<class ValueT>
class Vector_Abstract {

    /** Compute inner product with v
    @param v_s Other vector to compute inner product with.
        Indexing!
    @param invert If true, then compute 1/inner product. */
    virtual void apply(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,1> &out,
        FillType fill_type,
        bool invert = false) = 0;

};

/** Abstract "compressed" matrix, just enough to multiply by. */
template<class ValueT>
class Matrix_Abstract {

    BOOST_ENUM_VALUES(FillType, int,
        (zero_all) (0)    // Fill unused dimensions with 0
        (nan_all) (1)     // Fill unused dimensions with NaN
        (zero_some) (2)    // Zero out ONLY the dimensions we use
        (ignore) (3)  // Ignore unused imensions, don't touch them
    )



    /** Computes M*v_s, stores result in out_s. */
    virtual void apply(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type) = 0;

};
#endif


template<class ValueT>
class WeightedMatrix_Abstract {
    // TmpAlloc tmp;    // Sometimes, hold the things we're wrapping.

    virtual void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<ValueT,2> const &As,    // As(nvec, ndim)
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan,
        int invert=false) = 0;

public:
    /** Compute M * As */
    virtual void apply_M(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan,
        bool force_conservation=true) = 0;

    /** Compute wM * As */
    void apply_wM(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan)
    { apply_weight(0, As, out, fill_type, false); }

    /** Compute Mw * As */
    void apply_Mw(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan)
    { apply_weight(1, As, out, fill_type, false); }

    /** Compute 1. / (wM * As) */
    void apply_sM(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan)
    { apply_weight(0, As, out, fill_type, true); }

    /** Compute 1. / (Mw * As) */
    void apply_Ms(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan)
    { apply_weight(1, As, out, fill_type, true); }

};



};    // namespace
#endif    // guad
