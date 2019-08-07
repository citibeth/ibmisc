#ifndef SPSPARSE_TUPLELIST_HPP
#define SPSPARSE_TUPLELIST_HPP

#include <vector>
#include <array>
#include <ibmisc/netcdf.hpp>

namespace spsparse {

/** An N-dimensional generalization of Eigen::Triplet */
template<class IndexT, class ValT, int RANK>
class Tuple {
    std::array<IndexT, RANK> _index;
    ValT _value;
public:
    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & _index;
        ar & _value;
    }


    IndexT &index(int i)
        { return _index[i]; }
    IndexT const &index(int i) const
        { return _index[i]; }

    std::array<IndexT, RANK> &index()
        { return _index; }
    std::array<IndexT, RANK> const &index() const
        { return _index; }

    ValT &value()
        { return _value; }
    ValT const &value() const
        { return _value; }

    Tuple() {}    // WARNING: Uninitialized.  It's here for boost serialization

    Tuple(std::array<IndexT, RANK> const &index,
        ValT const &value)
        : _index(index), _value(value) {}

    // ----- For Eigen::SparseMatrix::setFromTriplets()
    ValT row() const
        { return index(0); }
    ValT col() const
        { return index(1); }

    bool operator<(Tuple<IndexT,ValT,RANK> const &other) const
        { return _index < other._index; }

    bool operator==(Tuple<IndexT,ValT,RANK> const &other) const
    {
        for (int i=0; i<RANK; ++i) if (_index[i] != other._index[i]) return false;
        return (_value == other._value);
    }

};

/** Serves as accumulator and iterable storage */
template<class IndexT, class ValT, int RANK>
class TupleList
{
public:
    // https://stackoverflow.com/questions/4353203/thou-shalt-not-inherit-from-stdvector
    typedef std::vector<Tuple<IndexT,ValT,RANK>> VectorT;
    VectorT tuples;

    // Stuff to make it an accumulator
    static const int rank = RANK;
    typedef IndexT index_type;
    typedef ValT val_type;
    
    std::array<long, rank> _shape;

public:
    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & tuples;
        ar & _shape;
    }

    // -----------------------------------------------------
    // https://stackoverflow.com/questions/7758580/writing-your-own-stl-container/7759622#7759622

    // These iterators work for Spsparse accumulators; but not for all STL stuff (eg std::sort).
    // For those, use TupleList.tuples.begin(), etc.
    struct iterator : public VectorT::iterator {
        static const int rank = RANK;
        typedef IndexT index_type;
        typedef ValT val_type;

        iterator() {}
        iterator(typename VectorT::iterator &&ii) : VectorT::iterator(std::move(ii)) {}
        iterator(const typename VectorT::iterator &ii) : VectorT::iterator(ii) {}

    };
    iterator begin()
        { return iterator(tuples.begin()); }
    iterator end()
        { return iterator(tuples.end()); }
    // -------------------------------------------------
    struct const_iterator : public VectorT::const_iterator {
        static const int rank = RANK;
        typedef IndexT index_type;
        typedef ValT val_type;

        const_iterator(typename VectorT::const_iterator &&ii) : VectorT::const_iterator(std::move(ii)) {}
    };
    const_iterator begin() const
        { return const_iterator(tuples.begin()); }
    const_iterator end() const
        { return const_iterator(tuples.end()); }

    // -------------------------------------------------

    TupleList()
        { _shape.fill(-1); }    // -1 means unlimited dimension

    TupleList(std::array<long,RANK> shape) : _shape(shape) {}

    // So this can serve as a Spsparse Accumulator
    void set_shape(std::array<long, rank> shape)
        { _shape = shape; }

    std::array<long,rank> const &shape() const
        { return _shape; }
    long shape(int i) const
        { return _shape[i]; }

    void add(std::array<index_type,rank> const &index, ValT const &value);

    // Forward methods to std::vector
    size_t size() const { return tuples.size(); }
    void clear() { tuples.clear(); }
    void reserve(size_t n) { tuples.reserve(n); }
    Tuple<IndexT,ValT,RANK> &operator[](int ix) { return tuples[ix]; }
    Tuple<IndexT,ValT,RANK> const &operator[](int ix) const { return tuples[ix]; }

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);
private:
    void nc_rw(netCDF::NcGroup *nc, char rw, std::string const &vname);
};

template<class IndexT, class ValT, int RANK>
void TupleList<IndexT,ValT,RANK>::add(std::array<index_type,rank> const &index, ValT const &value)
{
    // Check bounds
    for (int i=0; i<RANK; ++i) {
        if (_shape[i] >= 0 && (index[i] < 0 || index[i] >= _shape[i])) {
            std::ostringstream buf;
            buf << "Sparse index out of bounds: index=(";
            for (int j=0; j<RANK; ++j) {
                buf << index[j];
                buf << " ";
            }
            buf << ") vs. shape=(";
            for (int j=0; j<RANK; ++j) {
                buf << _shape[j];
                buf << " ";
            }
            buf << ")";
            (*ibmisc::ibmisc_error)(-1, buf.str().c_str());
        }
    }

    tuples.push_back(Tuple<IndexT,ValT,RANK>(index, value));
}

template<class IndexT, class ValT, int RANK>
void TupleList<IndexT,ValT,RANK>::nc_rw(
    netCDF::NcGroup *nc,
    char rw, std::string const &vname)
{
    netCDF::NcVar indices_v = nc->getVar(vname + ".indices");
    netCDF::NcVar vals_v = nc->getVar(vname + ".values");

    if (rw == 'w') {
        int const N = size();

        // Create in-memory data structure amenable to writing to disk quickly
        {std::vector<IndexT> indices;
            indices.reserve(N*RANK);
            for (auto &ii : tuples) {
                for (int k=0; k<RANK; ++k) indices.push_back(ii.index(k));
            }
            indices_v.putVar(&indices[0]);
        }

        // Write it out!
        {std::vector<ValT> vals;
            vals.reserve(N);
            for (auto &ii : tuples) {
                vals.push_back(ii.value());
            }
            vals_v.putVar(&vals[0]);    // Write to entire NetCDF variable directly from RAM
        }
    } else {    // rw == 'r'
        int const N = vals_v.getDim(0).getSize();

        // Create in-memory data structure amenable to writing to disk quickly
        std::vector<IndexT> indices;
        std::vector<ValT> vals;

        // Read it
        indices.resize(N*RANK);
        indices_v.getVar(&indices[0]);
        vals.resize(N);
        vals_v.getVar(&vals[0]);    // Write to entire NetCDF variable directly from RAM

        // Convert to TupleList
        // TODO: If we are more clever with iterators, we don't have to copy to
        //       a TupleList first.
        std::array<IndexT,RANK> ix;
        for (size_t i=0; i<vals.size(); ++i) {
            for (int k=0; k<RANK; ++k) ix[k] = indices[i*RANK+k];
            add(ix, vals[i]);
        }
        indices.clear();
        vals.clear();
    }
}

template<class IndexT, class ValT, int RANK>
void TupleList<IndexT,ValT,RANK>::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    std::vector<std::string> const dim_names();
    std::vector<netCDF::NcDim> dims;        // Dimensions in NetCDF

    dims = ibmisc::get_or_add_dims(ncio,
        {vname + ".nnz", vname + ".rank"},
        {size(),         RANK});    // length ign

    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    ibmisc::get_or_put_att(info_v, ncio.rw, "shape", "int64", _shape);

    get_or_add_var(ncio, vname + ".indices", ibmisc::get_nc_type<IndexT>(), dims);
    get_or_add_var(ncio, vname + ".values", ibmisc::get_nc_type<ValT>(), {dims[0]});
    ncio += std::bind(&TupleList<IndexT,ValT,RANK>::nc_rw, this, ncio.nc, ncio.rw, vname);

}



}   // namespace
#endif    // guard
