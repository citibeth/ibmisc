#include <ibmisc/indexing.hpp>

namespace ibmisc {

void Indexing::make_strides()
{
    data[indices[rank()-1]]._stride = 1;
    for (int d=rank()-2; d>=0; --d) {
        auto &data0(data[indices[d]]);
        auto &data1(data[indices[d+1]]);
        data0._stride = data1.stride() * data1.extent;
    }
}

Indexing::Indexing(
    std::vector<IndexingData> &&_data,
    std::vector<int> &&_indices)
: data(std::move(_data)),
    indices(std::move(_indices))
{ make_strides(); }

long Indexing::extent() const
{
    long ret = data[0].extent;
    for (int k=1; k<rank(); ++k) ret *= data[k].extent;
    return ret;
}

Indexing::Indexing(
    std::vector<std::string> const &_name,
    std::vector<long> const &_base,
    std::vector<long> const &_extent,
    std::vector<int> &&_indices) :
indices(std::move(_indices))
{
    std::vector<IndexingData> _data;
    for (size_t i=0; i<_base.size(); ++i)
        _data.push_back(IndexingData(_name[i], _base[i], _extent[i]));
    data = std::move(_data);
    make_strides();
}





void Indexing::ncio(
    NcIO &ncio,
    std::string const &vname)
{
    std::vector<std::string> name;
    std::vector<long> base;
    std::vector<long> extent;
    for (size_t i=0; i<rank(); ++i) {
        IndexingData const &data((*this)[i]);
        name.push_back(data.name);
        base.push_back(data.base);
        extent.push_back(data.extent);
    }

    auto info_v = get_or_add_var(ncio, vname, "int64", {});
    get_or_put_att(info_v, ncio.rw, "name", "string", name);
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
    std::vector<long> low;
    std::vector<long> high;
    for (size_t i=0; i<rank(); ++i) {
        DomainData const &data((*this)[i]);

        low.push_back(data.low);
        high.push_back(data.high);
    }

    auto info_v = get_or_add_var(ncio, vname, "int64", {});
    get_or_put_att(info_v, ncio.rw, "low", "int64", low);
    get_or_put_att(info_v, ncio.rw, "high", "int64", high);
}

bool in_domain(
    Domain const &domain,
    Indexing const &indexing,
    long ix)
{
    long tuple[domain.rank()];
    indexing.index_to_tuple(tuple, ix);
    return domain.in_domain(tuple);
}



};
