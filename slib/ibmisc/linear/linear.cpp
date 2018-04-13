#include <ibmisc/linear/linear.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <ibmisc/linear/compressed.hpp>

namespace ibmisc {
namespace linear {

void Weighted::to_coo(
    blitz::Array<int,1> &indices0,        // Must be pre-allocated(nnz)
    blitz::Array<int,1> &indices1,        // Must be pre-allocated(nnz)
    blitz::Array<double,1> &values) const
{
    long nnz = this->nnz();
    if (!indices0.data()) indices0.reference(blitz::Array<int,1>(nnz));
    if (!indices1.data()) indices1.reference(blitz::Array<int,1>(nnz));
    if (!values.data()) values.reference(blitz::Array<double,1>(nnz));

    if (indices0.extent(0) != nnz) (*ibmisc_error)(-1,
        "indices0[%d] must have dimension [%ld]",
        indices0.extent(0), nnz);

    if (indices1.extent(0) != nnz) (*ibmisc_error)(-1,
        "indices1[%d] must have dimension [%ld]",
        indices1.extent(0), nnz);

    if (values.extent(0) != nnz) (*ibmisc_error)(-1,
        "values[%d] must have dimension [%ld]",
        values.extent(0), nnz);

    _to_coo(indices0, indices1, values);
}


void Weighted::get_weights(
    int idim,    // 0=wM, 1=Mw
    blitz::Array<double,1> &w) const
{
    auto shape(this->shape());
    if (!w.data()) {
        w.reference(blitz::Array<double,1>(shape[idim]));
        w = 0;
    }

    if (w.extent(0) != shape[idim]) (*ibmisc_error)(-1,
        "weight%d[%d] must have dimnension [%d]",
        idim, w.extent(0), shape[idim]);

    _get_weights(idim, w);
}

void Weighted::ncio(NcIO &ncio, std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    if (ncio.rw == 'w') {
        LinearType _type = type;
        get_or_put_att_enum(info_v, ncio.rw, "type", _type);
    }
    get_or_put_att(info_v, ncio.rw, "conservative", conservative);
}

std::unique_ptr<Weighted> new_weighted(LinearType type)
{
    switch(type.index()) {
        case LinearType::EIGEN :
            return std::unique_ptr<Weighted>(new Weighted_Eigen);
        case LinearType::COMPRESSED :
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
