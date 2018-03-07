
// ==================================================================
class WeightedMatrix_Z : public WeightedMatrix_Abstract<double>
{
    std::array<ZSparseArray<int,double,1>, 2> weights;    // {wM, Mw}
    ZSparseArray<int,double,2> &M;

protected:
    /** Computes out = As * weights[dim] */
    void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<ValueT,2> const &As,    // As(nvec, ndim)
        blitz::Array<ValueT,2> &out,
        bool zero_out=true);

public:
    /** Computes out = M * As */
    void apply_M(
        blitz::Array<ValueT,2> const &As,
        blitz::Array<ValueT,2> &out,
        FillType fill_type=FillType::nan,
        bool force_conservation=true);

};





