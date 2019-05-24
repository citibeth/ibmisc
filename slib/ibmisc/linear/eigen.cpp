#include <ibmisc/linear/eigen.hpp>
#include <spsparse/eigen.hpp>
#include <ibmisc/netcdf.hpp>

using namespace spsparse;
using namespace blitz;

namespace ibmisc {
namespace linear {

// ----------------------------------------------------------------
static void mask_result(Weighted_Eigen::EigenDenseMatrixT &ret, blitz::Array<double,1> const &wB_b, double fill)
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
Weighted_Eigen::EigenDenseMatrixT Weighted_Eigen::apply_e(
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


    if (BvA.M->cols() != A.rows()) (*ibmisc_error)(-1,
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

void Weighted_Eigen::ncio(ibmisc::NcIO &ncio, std::string const &vname,
        std::array<std::string,2> dim_names)
{
    // Call to superclass
    Weighted::ncio(ncio, vname);

    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    //get_or_put_att(info_v, ncio.rw, "dim_names", "", dim_names);

    // If reading... overwrite user-supplied dim_names with what we found on disk.
    std::vector<std::string> data_v{dim_names[0], dim_names[1]};
    get_or_put_att(info_v, ncio.rw, "dim_names", "", data_v);
    if (ncio.rw == 'r') {
        dim_names[0] = data_v[0];
        dim_names[1] = data_v[1];
    }

    // Matches dimension name created by SparseSet.
    auto ncdims(ibmisc::get_or_add_dims(ncio,
        {dim_names[0] + ".dense_extent", dim_names[1] + ".dense_extent"},
        {dims[0]->dense_extent(), dims[1]->dense_extent()}));

    // --------- M
    if (ncio.rw == 'r') M.reset(new EigenSparseMatrixT);
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

void Weighted_Eigen::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    // Use standard dimnames for virtual function implementation.
    std::array<std::string,2> default_dim_names{vname+".dimB", vname+".dimA"};
    this->ncio(ncio, vname, default_dim_names);

    // Read/Write Dimensions, which might be shared with another Weighted_Eigen
    if (ncio.rw == 'r') {
        auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});

        // Retrieve actual dim_names this was written with.
        std::vector<std::string> data_v;
        get_or_put_att(info_v, ncio.rw, "dim_names", "", data_v);
        std::array<std::string,2> dim_names {data_v[0], data_v[1]};

        for (int i=0; i<2; ++i) {
            // Allocate dimension if not already done
            if (!dims[i]) dims[i] = tmp.newptr<SparseSetT>();

            // Read dimension if not already done
            if (dims[i]->dense_extent() == 0) dims[i]->ncio(ncio, dim_names[i]);
        }
    } else {
        for (int i=0; i<2; ++i) {
            // Write dimension if not already written
            if (ncio.nc->getVar(default_dim_names[i]).isNull()) {
                dims[i]->ncio(ncio, default_dim_names[i]);
            }
        }
    }
            
}
// ------------------------------------------------------
void Weighted_Eigen::apply_weight(
    int _dim,    // 0=B, 1=A
    blitz::Array<double,2> const &As,    // As(nvec, ndim)
    blitz::Array<double,1> &out,    // out(nvec)
    bool zero_out) const
{
    auto &weights(_dim == 0 ? wM : Mw);
    auto &dim(*dims[_dim]);
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));

    if (zero_out) out = 0;

    for (int j_d=0; j_d < dim.dense_extent(); ++j_d) {
        int const j_s = dim.to_sparse(j_d);
        for (int k=0; k<nvec; ++k) {
            out(k) += weights(j_d) * As(k,j_s);
        }
    }
}

/** Sparse shape of the matrix */
std::array<long,2> Weighted_Eigen::shape() const
    { return std::array<long,2>{dims[0]->sparse_extent(), dims[1]->sparse_extent()}; }


/** Compute M * As */
void Weighted_Eigen::apply_M(
    blitz::Array<double,2> const &A_s,
    blitz::Array<double,2> &B_s,
    AccumType accum_type,
    bool force_conservation) const
{
    // TODO: Re-do this method, to work without copying over the matrix.
    //       This would have to stop using Eigen's facilities

    // |j_s| = size of sparse input vector space (A_s)
    // |j_d] = size of dense input vector space (A_d)
    // |n| = number of variables being processed together

    // Allocate dense A matrix
    auto &bdim(*dims[0]);
    auto &adim(*dims[1]);
    int n_n = A_s.extent(0);

    // Densify the A matrix
    blitz::Array<double,2> A_d(n_n, adim.dense_extent());
    for (int j_d=0; j_d < adim.dense_extent(); ++j_d) {
        int const j_s = adim.to_sparse(j_d);
        for (int n=0; n < n_n; ++n) {
            A_d(n,j_d) = A_s(n,j_s);
        }
    }

    // Apply...
    auto B_d_eigen(apply_e(A_d, 0, force_conservation));    // Column major indexing

    switch(accum_type.value()) {
        case AccumType::REPLACE :
            for (int j_d=0; j_d < bdim.dense_extent(); ++j_d) {
                if (wM(j_d) == 0.) continue;    // Skip nullspace that crept into dense
                int j_s = bdim.to_sparse(j_d);
                for (int n=0; n < n_n; ++n) B_s(n,j_s) = B_d_eigen(j_d,n);
            }
        break;
        case AccumType::ACCUMULATE :
            for (int j_d=0; j_d < bdim.dense_extent(); ++j_d) {
                if (wM(j_d) == 0.) continue;    // Skip nullspace that crept into dense
                int j_s = bdim.to_sparse(j_d);
                for (int n=0; n < n_n; ++n) B_s(n,j_s) += B_d_eigen(j_d,n);
            }
        break;
        case AccumType::REPLACE_OR_ACCUMULATE :
            for (int j_d=0; j_d < bdim.dense_extent(); ++j_d) {
                if (wM(j_d) == 0.) continue;    // Skip nullspace that crept into dense
                int j_s = bdim.to_sparse(j_d);
                for (int n=0; n < n_n; ++n) {
                    auto &oval(B_s(n,j_s));
                    if (std::isnan(oval)) oval = B_d_eigen(j_d,n);
                    else oval += B_d_eigen(j_d,n);
                }
            }
        break;
    }
}

long Weighted_Eigen::nnz() const
{
    return M->nonZeros();
}

void Weighted_Eigen::_to_coo(
    blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
    blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
    blitz::Array<double,1> &values) const      // Must be pre-allocated(nnz)
{
    long nnz = this->nnz();
    long j = 0;
    for (auto ii(begin(*M)); ii != end(*M); ++ii) {
        indices0(j) = dims[0]->to_sparse(ii->index(0));
        indices1(j) = dims[1]->to_sparse(ii->index(1));
        values(j) = ii->value();
        ++j;
    }
}

void Weighted_Eigen::_get_weights(
    int idim,    // 0=wM, 1=Mw
    blitz::Array<double,1> &w) const
{
    auto &weights(idim == 0 ? wM : Mw);
    for (int j_d=0; j_d < dims[idim]->dense_extent(); ++j_d) {
        int const j_s = dims[idim]->to_sparse(j_d);
        w(j_s) += weights(j_d);
    }
}

}}    // namespace
