#include <ibmisc/linear/compressed.hpp>


namespace ibmisc {

static double const nan = std::numeric_limits<double>::quiet_NaN();

namespace linear {

void Weighted_Compressed::apply_weight(
    int dim,    // 0=B, 1=A
    blitz::Array<double,2> const &As,    // As(nvec, nA)
    blitz::Array<double,1> &out,          // out(nvec)
    bool zero_out)
{
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));

    if (zero_out) out = 0;

    for (auto ii(weights[dim].generator()); ++*ii; ) {
        for (int k=0; k<nvec; ++k) {
            out(k) += ii->value() * As(ii->index(0));
        }
    }
}

void Weighted_Compressed::apply_M(
    blitz::Array<double,2> const &As,    // As(nvec, nA)
    blitz::Array<double,2> &Bs,         // Bs(nvec, nB)
    FillType fill_type,
    bool force_conservation)
{
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));
    auto const nB(Bs.extent(1));

    // Zero out beforehand
    switch(fill_type.index()) {
        case FillType::zero_all :
            Bs = 0;
        break;
        case FillType::nan_all :
            Bs = nan;
        break;
        case FillType::zero_some :
            for (auto ii(weights[0].generator()); ++*ii; ) {
                for (int k=0; k<nvec; ++k) Bs(k,ii->index(0)) = 0;
            }
        break;
    }

    // Multiply
    for (auto ii(M.generator()); ++*ii; ) {
        auto const i(ii->index(0));

        for (int k=0; k<nvec; ++k) {
            auto &B_ki(Bs(k,i));
            auto val(ii->value() * As(ii->index(1)));

            B_ki = ((fill_type == FillType::nan_all) && std::isnan(B_ki)
                ? val : B_ki + val);
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
        for (auto ii(weights[0].generator()); ++*ii; ) {
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

}}    // namespace
