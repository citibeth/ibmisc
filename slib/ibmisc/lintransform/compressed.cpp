void Compressed::apply_weight(
    int dim,    // 0=B, 1=A
    blitz::Array<ValueT,2> const &As,    // As(nvec, nA)
    blitz::Array<ValueT,2> &out,          // out(nvec)
    bool zero_out=true)
{
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));

    if (zero_out) out = 0;

    for (auto ii((dim == 0 ? wM : Mw).generator()); ++*ii; ) {
        for (int k=0; k<nvec; ++k) {
            out(k) += ii->value() * As(ii->index(0));
        }
    }
}

void Compressed::apply_M(
    blitz::Array<double,2> const &As,    // As(nvec, nA)
    blitz::Array<double,2> &out,         // out(nvec, nB)
    FillType fill_type=FillType::nan,
    bool correct_conservation=false)
{
    auto const nvec(As.extent(0));
    auto const nA(As.extent(1));
    aout const nB(Bs.extent(1));

    // Zero out beforehand
    switch(fill_type.index()) {
        case FillType::zero_all :
            out = 0;
        break;
        case FillType::nan_all :
            out = nan;
        break;
        case FillType::zero_some :
            for (auto ii(weights[0].generator()); ++*ii; ) {
                for (int j=0; j<nvec; ++j) out(nvec,ii->index(0)) = 0;
            }
        break;
    }

    // Multiply
    for (auto ii(M.generator()); ++*ii; ) {
        auto const i(ii->index(0));

        for (int k=0; k<nvec; ++k) {
            auto &B_ki(out(k,i));
            auto val(ii->value() * As(ii->index(1)));

            B_ki = ((fill_type == FillType::nan_all) && std::isnan(B_ki)
                ? val : B_ki + val);
        }
    }

    // correct_conservation = (force_conservation && !conservative)
    if (correct_conservation) {
        // Compute correction factor for each variable
        blitz::Array<double,1> wA(nvec);
        blitz::Array<double,1> wB(nvec);

        apply_weight(0, out, wB, true);
        apply_weight(1, As, wA, true);

        blitz::Array<double,1> factor(nvec);
        for (int k=0; k<nvec; ++k) factor(k)=wA(k)/wB(k);

        // Multiply by correction factor
        for (auto ii(wM.generator()); ++*ii; ) {
            auto const i(ii->index[0]);
            for (int k=0; k<nvec; ++k) {
                out(k,i) *= factor(k);
            }
        }
    }
}

void Compressed::ncio(NcIO &ncio, std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att(info_v, ncio.rw, "conservative", conservative);

    weights[0].ncio(ncio, vname + ".wM");
    M.ncio(ncio, vname + ".M");
    weights[1].ncio(ncio, vname + ".Mw");
}
