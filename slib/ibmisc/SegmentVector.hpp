#include <vector>

namespace ibmisc {

/** Stores a series of sub-vectors, logically combining them into one vector.
Unit test is in test_iter.cpp*/
template<class ItemT>
class SegmentIterator
{
    typename std::vector<std::vector<ItemT>>::iterator const _end0;
    typename std::vector<std::vector<ItemT>>::iterator _i0;
    typename std::vector<ItemT>::iterator _i1;

    void normalize()
    {
        if (_i0 == _end0) return;
        while (_i1 == _i0->end()) {
            ++_i0;
            _i1 = _i0->begin();     // Safe because of sentinel
            if (_i0 == _end0) break;
        }
    }

public:

    SegmentIterator(
        typename std::vector<std::vector<ItemT>>::iterator end0,
        typename std::vector<std::vector<ItemT>>::iterator i0,
        typename std::vector<ItemT>::iterator i1)
    : _end0(end0), _i0(i0), _i1(i1)
    {
        normalize();        // Seek to first valid value
    }

    ItemT &operator*()
        { return *_i1; }

    SegmentIterator& operator++()
    {
        ++_i1;
        normalize();
    }

    bool operator==(SegmentIterator const &rhs) const
        { return (_i0 == rhs._i0) && (_i1 == rhs._i1); }
    bool operator!=(const SegmentIterator& rhs) const
        {return !this->operator==(rhs); }
};

template<class ItemT>
class SegmentVector {
    std::vector<std::vector<ItemT>> segments;

public:

    // ------------------------------------------------------------------------
    typedef typename std::vector<std::vector<ItemT>>::iterator segment_iterator;
    typedef typename std::vector<std::vector<ItemT>>::const_iterator segment_const_iterator;

    segment_iterator segment_begin() { return segments.begin(); }
    segment_iterator segment_end() { return std::prev(segments.end()); }    // Sentinel
    segment_const_iterator segment_begin() const { return segments.begin(); }
    segment_const_iterator segment_end() const { return std::prev(segments.end()); }

    std::vector<ItemT> &segment(int ix) { return segments[ix]; }
    // ------------------------------------------------------------------------

    typedef SegmentIterator<ItemT> iterator;
    typedef SegmentIterator<const ItemT> const_iterator;

    SegmentVector(size_t _nsegments) : segments(_nsegments+1) {}    // sentinel

    size_t nsegments() { return segments.size()-1; }        // sentinel

    SegmentIterator<ItemT> begin()
        { return SegmentIterator<ItemT>(segments.end(), segment_begin(), segment_begin()->begin()); }
    SegmentIterator<ItemT> end()
        { return SegmentIterator<ItemT>(segments.end(), segments.end(), segment_end()->begin()); }


    // -------- These could be optimized...
    size_t size()
    {
        size_t ret = 0;
        for (auto i0=segments.begin(); i0 != segments.end(); ++i0) ret += i0->size();
        return ret;
    }

    ItemT &operator[](int ix)
    {
        // Figure out which segment it's in
        size_t iseg=0;
        size_t seg_start=0;
        for (; iseg < nsegments();
            seg_start += segment(iseg).size(), ++iseg)
        {
            size_t const local_ix = ix-seg_start;
            if (local_ix < segment(iseg).size()) return segment(iseg)[local_ix];
        }

        // Out of range; let underlying STL throw the exception
        return segment(nsegments())[0];
    }

};

}
