#include <spsparse/eigen.hpp>
#include <spsparse/blitz.hpp>
#include <spsparse/SparseSet.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <ibmisc/linear/eigen.hpp>


namespace ibmisc {

static double const nan = std::numeric_limits<double>::quiet_NaN();

namespace linear {

void Weighted_Compressed::set_shape(std::array<long,2> _shape)
{
    M.accum().set_shape(_shape);
}

std::array<long,2> Weighted_Compressed::shape() const
{
    return M.shape();
}


void Weighted_Compressed::apply_weight(
    int dim,    // 0=B, 1=A
    blitz::Array<double,2> const &As,    // As(nvec, nA)
    blitz::Array<double,1> &out,          // out(nvec)
    bool zero_out) const
{
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));

    if (zero_out) out = 0;

    for (auto ii(weights[dim].generator()); ++ii; ) {
        for (int k=0; k<nvec; ++k) {
            out(k) += ii->value() * As(k,ii->index(0));
        }
    }
}


void Weighted_Compressed::apply_M(
    blitz::Array<double,2> const &As,    // As(nvec, nA)
    blitz::Array<double,2> &Bs,         // Bs(nvec, nB)
    AccumType accum_type,
    bool force_conservation) const
{
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));
    auto const nB(Bs.extent(1));

    // Prepare ActiveSpace, based on accum_type
    switch(accum_type.index()) {
        case AccumType::REPLACE :
            for (auto ii(weights[0].generator()); ++ii; ) {
                for (int k=0; k<nvec; ++k) Bs(k,ii->index(0)) = 0;
            }
        break;
        case AccumType::REPLACE_OR_ACCUMULATE :
            for (auto ii(weights[0].generator()); ++ii; ) {
                for (int k=0; k<nvec; ++k) {
                    auto &Bs_ki(Bs(k,ii->index(0)));
                    if (std::isnan(Bs_ki)) Bs_ki = 0;
                }
            }
        break;
    }

    // Multiply
    for (auto ii(M.generator()); ++ii; ) {
        auto const i(ii->index(0));

        for (int k=0; k<nvec; ++k) {
            Bs(k,i) += ii->value() * As(k,ii->index(1));
        }
    }

    if (force_conservation && !conservative) {
        // Compute correction factor for each variable
        blitz::Array<double,1> wA(nvec);
        blitz::Array<double,1> wB(nvec);

        apply_weight(0, Bs, wB, true);
        apply_weight(1, As, wA, true);

        blitz::Array<double,1> factor(nvec);
        for (int k=0; k<nvec; ++k) factor(k)=wA(k)/wB(k);

        // Multiply by correction factor
        for (auto ii(weights[0].generator()); ++ii; ) {
            auto const i(ii->index(0));
            for (int k=0; k<nvec; ++k) {
                Bs(k,i) *= factor(k);
            }
        }
    }
}

void Weighted_Compressed::ncio(NcIO &ncio, std::string const &vname)
{
    // Call to superclass
    Weighted::ncio(ncio, vname);

    weights[0].ncio(ncio, vname + ".wM");
    M.ncio(ncio, vname + ".M");
    weights[1].ncio(ncio, vname + ".Mw");
}

// ======================================================
Weighted_Compressed compress(Weighted_Eigen &eigen)
{
    Weighted_Compressed ret(eigen.scaled);
    ret.conservative = eigen.conservative;

    spsparse::spcopy(
        spsparse::accum::to_sparse(
            std::array<Weighted_Eigen::SparseSetT *,1>{eigen.dims[0]},
        ret.weights[0].accum()),
        eigen.wM);

    spsparse::spcopy(
        spsparse::accum::to_sparse(eigen.dims,
        ret.M.accum()),
        *eigen.M);

    spsparse::spcopy(
        spsparse::accum::to_sparse(
            std::array<Weighted_Eigen::SparseSetT *,1>{eigen.dims[1]},
        ret.weights[1].accum()),
        eigen.Mw);

    return ret;
}

long Weighted_Compressed::nnz() const
    { return M.nnz(); }

void Weighted_Compressed::_to_coo(
    blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
    blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
    blitz::Array<double,1> &values) const      // Must bepre-allocated(nnz)
{
    long nnz = this->nnz();
    long j = 0;
    for (auto ii(M.generator()); ++ii; ) {
        indices0(j) = ii->index(0);
        indices1(j) = ii->index(1);
        values(j) = ii->value();
        ++j;
    }
}

void Weighted_Compressed::_get_weights(
    int idim,    // 0=wM, 1=Mw
    blitz::Array<double,1> &w) const
{
    for (auto ii(weights[idim].generator()); ++ii; ) {
        w(ii->index(0)) += ii->value();
    }
}


}}    // namespace
