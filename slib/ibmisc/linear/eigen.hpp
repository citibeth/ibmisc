#ifndef IBMISC_LINEAR_EIGEN_HPP
#define IBMISC_LINEAR_EIGEN_HPP

#include <memory>
#include <Eigen/SparseCore>
#include <spsparse/SparseSet.hpp>
#include <ibmisc/linear/linear.hpp>

namespace ibmisc {
namespace linear {

/** Return value of a sparse matrix */
struct Weighted_Eigen : public Weighted {
    typedef spsparse::SparseSet<long,int> SparseSetT;
    typedef Eigen::SparseMatrix<double,0,int> EigenSparseMatrixT;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> EigenColVectorT;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic> EigenRowVectorT;

    ibmisc::TmpAlloc tmp;            // Stuff that needs same lifetime as Eigen

    std::array<std::string,2> dim_names;
    std::array<SparseSetT *,2> dims;            // Dense-to-sparse mapping for the dimensions.  (Store the originals in tmp if you need to)

    // If M=BvA, then wM = wBvA = area of B cells
    blitz::Array<double,1> wM;           // Dense indexing


    std::unique_ptr<EigenSparseMatrixT> M;    // Dense indexing

    // Area of A cells
    blitz::Array<double,1> Mw;

    /** Construct with shared dimensions from elsewhere */
    Weighted_Eigen(std::array<std::string,2> const &_dim_names, std::array<SparseSetT *,2> _dims, bool conservative=true)
        : Weighted(LinearType::EIGEN, conservative), dim_names(_dim_names), dims(_dims) {}

    /** Construct with internal dimensions.
    @param dim_names Starting with a dot means, add vname in ncio().
    NOTE: dim_names might be changed later, as we read from disk. */
    Weighted_Eigen(std::array<std::string,2> const &_dim_names={".dimB",".dimA"}, bool conservative=true) : Weighted(LinearType::EIGEN, conservative), dim_names(_dim_names)
    {
        dims[0] = tmp.newptr<SparseSetT>();
        dims[1] = tmp.newptr<SparseSetT>();
    }

    /** Applies a regrid matrix.
    Nominally computes B{in} = BvA{ij} * A{jn}
    (where BvA is this)
    If the matrix ix non-conservative, it also accounts for that.

        |i| = Size of B vector space
        |j| = Size of A vector space
        |n| = Number of vectors being transformed

    @param A The values to regrid, as a series of Eigen column vectors.
    @return Eigen type
    */
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> apply_e(
        // Eigen const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable (_b means dense here)
        double fill = std::numeric_limits<double>::quiet_NaN(),    // Fill value for cells not in BvA matrix
        bool force_conservation=true) const;     // Set if you want apply_e() to conserve, even if !M->conservative


    /** Apply to multiple variables
    @return Blitz type */
    blitz::Array<double,2> apply(
        // Eigen const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
        double fill,    // Fill value for cells not in BvA matrix
        bool force_conservation,
        ibmisc::TmpAlloc &tmp) const;


    /** Apply to a single variable */
    blitz::Array<double,1> apply(
        // Eigen const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,1> const &A_b,       // A_b{j} One variable
        double fill,    // Fill value for cells not in BvA matrix
        bool force_conservation,
        ibmisc::TmpAlloc &tmp) const;

    /** NOTE: this->dims must already be allocated and read, if you are reading. */
    void ncio_nodim(ibmisc::NcIO &ncio, std::string const &vname);

    /** Read/write to NetCDF */
    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    // =============== Implement linear::Weighted
    void apply_weight(
        int dim,    // 0=B, 1=A
        blitz::Array<double,2> const &As,    // As(nvec, ndim)
        blitz::Array<double,1> &out,
        bool zero_out=true) const;

    /** Sparse shape of the matrix */
    std::array<long,2> shape() const;

    /** Compute M * As */
    void apply_M(
        blitz::Array<double,2> const &As,
        blitz::Array<double,2> &out,
        AccumType accum_type=AccumType::REPLACE,
        bool force_conservation=true) const;

};


}}    // namespace
#endif
