#include <ibmisc/lintransform/eigen.hpp>

namespace ibmisc {
namespace lintransform {

// ----------------------------------------------------------------
static void mask_result(EigenDenseMatrixT &ret, blitz::Array<double,1> const &wB_b, double fill)
{
    int nB = ret.rows();    // == wB_b.extent(0)
    int nvar = ret.cols();

    // Mask out cells that slipped into the output because they were
    // in the SparseSet; but don't actually get any contribution.
    for (int i=0; i<nB; ++i) {
        if (wB_b(i) != 0.) continue;
        for (int n=0; n<nvar; ++n) ret(i,n) = fill;
    }

}
// -----------------------------------------------------------------------
EigenDenseMatrixT Weighted_Eigen::apply_e(
    // Eigen const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
    double fill,     // Fill value for cells not in BvA matrix
    bool force_conservation) const
{
    auto &BvA(*this);

    // A{jn}   One col per variable
    int nvar = A_b.extent(0);
    int nA = A_b.extent(1);

    Eigen::Map<EigenDenseMatrixT> const A(
        const_cast<double *>(A_b.data()), nA, nvar);

    // |i| = size of output vector space (B)
    // |j| = size of input vector space (A)
    // |n| = number of variables being processed together


    if (BvA.M->cols() != A.rows()) (*icebin_error)(-1,
        "BvA.cols=%d does not match A.rows=%d", BvA.M->cols(), A.rows());

    // Apply initial regridding.
    EigenDenseMatrixT B0(*BvA.M * A);        // B0{in}

    // Only apply conservation correction if all of:
    //   a) Matrix is smoothed, so it needs a conservation correction
    //   b) User requested conservation be maintained
    if (BvA.conservative || !force_conservation) {
        // Remove cells not in the sparse matrix
        mask_result(B0, BvA.wM, fill);
        return B0;
    }
    // -------------- Apply the Conservation Correction
    // Integrate each variable of input (A) over full domain
    auto &wA_b(BvA.Mw);
    Eigen::Map<EigenRowVectorT> const wA(const_cast<double *>(wA_b.data()), 1, wA_b.extent(0));
    typedef Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> EigenArrayT;
    EigenArrayT TA((wA * A).array());        // TA{n} row array

    // Integrate each variable of output (B) over full domain
    auto &wB_b(BvA.wM);
    int nB = wB_b.extent(0);
    Eigen::Map<EigenRowVectorT> const wB(const_cast<double *>(wB_b.data()), 1, wB_b.extent(0));
    EigenArrayT TB((wB * B0).array());
    EigenArrayT TB_inv(1. / TB);    // TB_inv{n} row array

    // Factor{nn}: Conservation correction for each variable.

    auto Factor(TA * TB_inv);    // Factor{n} row array

    std::cout << "-------- Eigen::apply() conservation" << std::endl;
    std::cout << "    |input|    = " << TA << std::endl;
    std::cout << "    |output|   = " << TB << std::endl;
    std::cout << "    correction = " << Factor << std::endl;

    EigenDenseMatrixT ret(B0 * Factor.matrix().asDiagonal());    // ret{in}
    // Remove cells not in the sparse matrix
    mask_result(B0, BvA.wM, fill);

    return ret;
}

blitz::Array<double,2> Weighted_Eigen::apply(
    // Eigen const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
    double fill,    // Fill value for cells not in BvA matrix
    bool force_conservation,
    ibmisc::TmpAlloc &tmp) const
{
    return spsparse::to_blitz<double>(apply_e(A_b, fill), tmp);
}


/** Apply to a single variable */
blitz::Array<double,1> Weighted_Eigen::apply(
    // Eigen const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,1> const &A_b,       // A_b{j} One variable
    double fill,    // Fill value for cells not in BvA matrix
    bool force_conservation,
    ibmisc::TmpAlloc &tmp) const
{
    auto A_b2(ibmisc::reshape<double,1,2>(A_b, {1, A_b.shape()[0]}));
    auto ret2(spsparse::to_blitz(apply_e(A_b2, fill), tmp));
    return ibmisc::reshape<double,2,1>(ret2, {ret2.shape()[1]});
}
// -----------------------------------------------------------------------
// ---------------------------------------------------------


void Weighted_Eigen::ncio(ibmisc::NcIO &ncio,
    std::string const &vname,
    std::array<std::string,2> dim_names)
{
    // Matches dimension name created by SparseSet.
    auto ncdims(ibmisc::get_or_add_dims(ncio,
        {dim_names[0] + ".dense_extent", dim_names[1] + ".dense_extent"},
        {dims[0]->dense_extent(), dims[1]->dense_extent()}));

    // --------- M
    ncio_eigen(ncio, *M, vname + ".M");

    // ---- Mw
    ncio_blitz_alloc<double,1>(ncio, Mw, vname + ".Mw", get_nc_type<double>(),
        {ncdims[1]});

    // ----------- wM
    std::string matrix_name(dim_names[0] + "v" + dim_names[1]);
    ncio_blitz_alloc<double,1>(ncio, wM, vname+".wM", get_nc_type<double>(),
        {ncdims[0]});

    netCDF::NcVar ncvar = ncio.nc->getVar(vname + ".M.info");
    get_or_put_att(ncvar, ncio.rw,
        "conservative", get_nc_type<bool>(), &conservative, 1);

}

// ------------------------------------------------------
void apply_weight(
    int dim,    // 0=B, 1=A
    blitz::Array<ValueT,2> const &As,    // As(nvec, ndim)
    blitz::Array<ValueT,2> &out,
    FillType fill_type=FillType::nan,
    int invert=false)
{
    TODO... invent new code
    REMEMBER, this must work with SPARSE indexing!
}

/** Sparse shape of the matrix */
std::array<long,2> shape()
{ TODO }


/** Compute M * As */
void apply_M(
    blitz::Array<ValueT,2> const &As,
    blitz::Array<ValueT,2> &out,
    FillType fill_type=FillType::nan,
    bool force_conservation=true)
{
    TODO... move code from icebin_cython.cpp
    REMEMBER, this must work with SPARSE indexing!
}




}}    // namespace
