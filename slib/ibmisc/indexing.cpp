#include <ibmisc/indexing.hpp>

namespace ibmisc {

void IndexingBase::make_strides()
{
    data[indices[rank()-1]].stride = 1;
    for (int d=rank()-2; d>=0; --d) {
        auto &data0(data[indices[d]]);
        auto &data1(data[indices[d+1]]);
        data0.stride = data1.stride * data1.extent;
    }
}

IndexingBase::IndexingBase(
    std::vector<IndexingData> &&_data,
    std::vector<int> &&_indices)
: data(std::move(_data)),
    indices(std::move(_indices))
{ make_strides(); }

long IndexingBase::extent() const
{
    long ret = extent[0];
    for (int k=1; k<rank(); ++k) ret *= extent[k];
    return ret;
}



void IndexingBase::ncio(
    NcIO &ncio,
    std::string const &vname)
{
    std::vector<std::string> name;
    std::vector<long> base;
    std::vector<long> extent;
    for (size_t i=0; i<rank(); ++i) {
        IndexingData &data((*this)[i]);
        name.push_back(data.name);
        base.push_back(data.base);
        extent.push_back(data.extent);
    }

    auto info_v = get_or_add_var(ncio, vname, "int64", {});
    get_or_put_att(info_v, ncio.rw, "name", netCDF::ncString, name);
    get_or_put_att(info_v, ncio.rw, "base", "int64", base);
    get_or_put_att(info_v, ncio.rw, "extent", "int64", extent);
    get_or_put_att(info_v, ncio.rw, "indices", "int64", indices);
    make_strides();
}

// ------------------------------------------------

void Domain::ncio(
    NcIO &ncio,
    std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname, "int64", {});
    get_or_put_att(info_v, ncio.rw, "low", "int64", low);
    get_or_put_att(info_v, ncio.rw, "high", "int64", high);
}

bool in_domain(
    DomainBase const &domain,
    BaseIndexing const &indexing,
    IndexT ix)
{
    long tuple[domain.rank()];
    indexing.index_to_tuple(tuple, ix);
    return domain.in_domain(tuple);
}



};
