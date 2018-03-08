#include <ibmisc/linear/linear.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <ibmisc/linear/compressed.hpp>

namespace ibmisc {
namespace linear {

void Weighted::ncio(NcIO &ncio, std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att(info_v, ncio.rw, "conservative", conservative);
}

std::unique_ptr<Weighted> new_weighted(LinearType type)
{
    switch(type.index()) {
        case Linear::EIGEN :
            return std::unique_ptr<Weighted>(new Weighted_Eigen);
        case Linear::COMPRESSED :
            return std::unique_ptr<Weighted>(new Weighted_Compressed);
        default:
            (*ibmisc_error)(-1,
                "Unrecognized LinearType = %d", type.index());
    }
}

std::unique_ptr<Weighted> nc_read_weighted(netCDF::NcGroup *nc, std::string const &vname)
{
    NcIO ncio(nc, 'r');    // Dummy

    std::string vn(vname + ".info");
    auto info_v = get_or_add_var(ncio, vn, "int", {});

    LinearType type;
    get_or_put_att_enum(info_v, ncio.rw, "type", type);
    auto ret(new_weighted(type));
    ret->ncio(ncio, vname);
    return ret;
}


}}    // namespace
