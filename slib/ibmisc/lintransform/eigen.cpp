#include <ibmisc/lintransform/eigen.hpp>

namespace ibmisc {
namespace lintransform {

static std::unique_ptr<Eigen> compute_AEvI(
    IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params,
    blitz::Array<double,1> const *elevmaskI,
    UrAE const &AE)
{
printf("BEGIN compute_AEvI scale=%d correctA=%d\n", params.scale, params.correctA);
    std::unique_ptr<Eigen> ret(new Eigen(dims, true));
    SparseSetT * const dimA(ret->dims[0]);
    SparseSetT * const dimI(ret->dims[1]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(AE.nfull);
    if (dimI) dimI->set_sparse_extent(regridder->nI());
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvI_m(
        // Only includes ice model grid cells with ice in them.
        std::bind(&IceRegridder::GvI, regridder, _1, elevmaskI),
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimI}, '.');
    MakeDenseEigenT ApvG_m(        // _m ==> type MakeDenseEigenT
        // Only includes ice model grid cells with ice in them.
        AE.GvAp,
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimA}, 'T');

    // ----- Convert to Eigen and multiply
    auto ApvG(ApvG_m.to_eigen());
    auto GvI(GvI_m.to_eigen());
    auto sGvI(sum_to_diagonal(GvI, 0, '-'));

    std::unique_ptr<EigenSparseMatrixT> ApvI(
        new EigenSparseMatrixT(ApvG * sGvI * GvI));
    ret->Mw.reference(sum(*ApvI, 1, '+'));    // Area of I cells

    // ----- Apply final scaling, and convert back to sparse dimension
    if (params.correctA) {
        // ----- Compute the final weight matrix
        auto wAvAp(MakeDenseEigenT(                   // diagonal
            AE.sApvA,
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
            {dimA, dimA}, '.').to_eigen());

        auto wApvI(sum_to_diagonal(*ApvI, 0, '+'));        // diagonal

        EigenSparseMatrixT wAvI(wAvAp * wApvI);    // diagonal...

        // +correctA: Weight matrix in A space
        ret->wM.reference(sum(wAvI, 0, '+'));    // Area of A cells

        // Compute the main matrix
        auto sAvAp(sum_to_diagonal(wAvAp, 0, '-'));
        if (params.scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI(sum_to_diagonal(wApvI, 0, '-'));

            ret->M.reset(new EigenSparseMatrixT(
                sAvAp * sApvI * *ApvI));    // AvI_scaled
        } else {
            // Should be like this for test_conserv.py
            // Note that sAvAp * sApvI = [size (weight) of grid cells in A]
            ret->M = std::move(ApvI);
        }

    } else {

        // ----- Compute the final weight matrix
        // ~correctA: Weight matrix in Ap space
        auto wApvI_b(sum(*ApvI, 0, '+'));
        ret->wM.reference(wApvI_b);

        if (params.scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI(diag_matrix(wApvI_b, '-'));

            ret->M.reset(new EigenSparseMatrixT(
                sApvI * *ApvI));    // ApvI_scaled
        } else {
            ret->M = std::move(ApvI);
        }
    }

printf("END compute_AEvI\n");
    return ret;
}
// ---------------------------------------------------------


void Eigen::ncio(ibmisc::NcIO &ncio,
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

}}    // namespace
