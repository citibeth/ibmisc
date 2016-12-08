#include <ibmisc/indexing.hpp>
#include <ibmisc/memory.hpp>

using namespace std;

namespace ibmisc {

bool IndexingData::operator==(IndexingData const &other) {
    return (name == other.name && base == other.base && extent == other.extent);
}

void Indexing::make_strides()
{
    if (rank() == 0) return;

    data[_indices[rank()-1]]._stride = 1;
    for (int d=rank()-2; d>=0; --d) {
        auto &data0(data[_indices[d]]);
        auto &data1(data[_indices[d+1]]);
        data0._stride = data1.stride() * data1.extent;
    }
}

Indexing::Indexing(
    std::vector<IndexingData> &&_data,
    std::vector<int> &&indices)
: data(std::move(_data)),
    _indices(std::move(indices))
{ make_strides(); }

bool Indexing::operator==(Indexing const &other) {
    if (data.size() != other.data.size()) return false;
    if (_indices != other._indices) return false;
    for (size_t i=0; i<data.size(); ++i) {
        if (!(data[i] == other.data[i])) return false;
    }
    return true;
}


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
    std::vector<int> &&indices) :
_indices(std::move(indices))
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
    auto &name(ncio.tmp.make<vector<std::string>>());
    auto &base(ncio.tmp.make<vector<long>>());
    auto &extent(ncio.tmp.make<vector<long>>());

    if (ncio.rw == 'w') for (size_t i=0; i<rank(); ++i) {
        IndexingData const &data((*this)[i]);
        name.push_back(data.name);
        base.push_back(data.base);
        extent.push_back(data.extent);
    }

    auto info_v = get_or_add_var(ncio, vname, "int64", {});
    get_or_put_att(info_v, ncio.rw, "name", "string", name);
    get_or_put_att(info_v, ncio.rw, "base", "int64", base);
    get_or_put_att(info_v, ncio.rw, "extent", "int64", extent);
    get_or_put_att(info_v, ncio.rw, "indices", "int64", _indices);


    if (ncio.rw == 'r') {
        data.clear();
        for (size_t i=0; i<name.size(); ++i) {
            data.push_back(IndexingData(name[i], base[i], extent[i]));
        }
        make_strides();
        ncio.tmp.free();    // Free low and high if reading
    }
}
// ------------------------------------------------

bool DomainData::operator==(DomainData const &other) {
    return (low == other.low && high == other.high);
}
bool Domain::operator==(Domain const &other) {
    if (data.size() != other.data.size()) return false;
    for (size_t i=0; i<data.size(); ++i) {
        if (!(data[i] == other.data[i])) return false;
    }
    return true;
}

void Domain::ncio(
    NcIO &ncio,
    std::string const &vname)
{
    auto &low(ncio.tmp.make<vector<long>>());
    auto &high(ncio.tmp.make<vector<long>>());

    if (ncio.rw == 'w') for (size_t i=0; i<rank(); ++i) {
        DomainData const &data((*this)[i]);
        low.push_back(data.low);
        high.push_back(data.high);
    }

    auto info_v = get_or_add_var(ncio, vname, "int64", {});
    get_or_put_att(info_v, ncio.rw, "low", "int64", low);
    get_or_put_att(info_v, ncio.rw, "high", "int64", high);

    if (ncio.rw == 'r') {
        data.clear();
        for (size_t i=0; i<low.size(); ++i) {
            data.push_back(DomainData(low[i], high[i]));
        }
        ncio.tmp.free();    // Free low and high if reading
    }
}

bool in_domain(
    Domain const *domain,
    Indexing const *indexing,
    long ix)
{
    long tuple[domain->rank()];
    indexing->index_to_tuple(tuple, ix);
    return domain->in_domain(tuple);
}



};
