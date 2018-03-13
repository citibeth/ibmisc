
namespace ibmisc {
namespace cython {


PyObject *linear_Weighted_shape(
    linear::Weighted *M)
{
    return py_build_tuple<pytype_long, long, 2>(M->shape());
}

PyObject *linear_Weighted_apply_weight(
    linear::Weighted *self,
    int dim,
    PyObject *A_s_py)            // A_b{nj_s} One row per variable
{
    auto shape(self->shape());

    // |j_s| = size of sparse input vector space (A_s)
    // |j_d] = size of dense input vector space (A_d)
    // |n| = number of variables being processed together

    // Convert input arrays
    auto A_s(np_to_blitz<double,2>(A_s_py, "A", {-1,-1}));    // Sparse indexed dense vector
    int n_n = A_s.extent(0);

    // Allocate output
    PyObject *B_s_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,2>{n_n});
    auto B_s(np_to_blitz<double1>(B_s_py, "B_s_py", {-1,-1}));
    B_s = nan;

    // Run it!
    self->apply_weight(dim, A_s, B_s, true);

    return B_s_py;
}

PyObject *linear_Weighted_apply_M(
    linear::Weighted *self,
    PyObject *A_s_py,            // A_b{nj_s} One row per variable
    double fill,
    // std::string const &saccum_type,
    bool force_conservation)
{
    // Convert input arrays
    auto A_s(np_to_blitz<double,2>(A_s_py, "A", {-1,-1}));    // Sparse indexed dense vector
    int n_n = A_s.extent(0);

    // Convert accum_type
    // auto accum_type(parse_enum<AccumType>(saccum_type));

    // Allocate output
    PyObject *B_s_py = ibmisc::cython::new_pyarray<double,2>(
        std::array<int,2>{n_n, shape[0]});
    auto B_s(np_to_blitz<double,2>(B_s_py, "B_s_py", {-1,-1}));
    B_s = fill

    // Run it!
    self->apply_M(A_s, B_s, AccumType::REPLACE, force_conservation);

    return B_s_py;
}


// -----------------------------------------
typedef spsparse::MakeDenseEigen<sparse_index_type, val_type, 0, dense_index_type> MakeDenseEigenT;
template<int RANK>
    using TupleListT = MakeDenseEigenT::TupleListT<RANK>;
template<int RANK>
    using DenseArrayT = blitz::Array<val_type,RANK>;
typedef MakeDenseEigenT::SparseSetT SparseSetT;
typedef MakeDenseEigenT::EigenSparseMatrixT EigenSparseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic> EigenRowVectorT;
typedef Eigen::DiagonalMatrix<val_type, Eigen::Dynamic> EigenDiagonalMatrixT;
// -----------------------------------------

void sample_makedense(MakeDenseEigenT::AccumT &&accum)
{

    accum.add({0,0},.5);
    accum.add({10,10},.5);
    accum.add({20,20},.5);


    //accum.add({30,30},.5);    // Nullspace
    accum.add({00,40},.5);
    accum.add({10,40},.5);
    accum.add({20,40},.5);
    //accum.add({30,40},.5);
}


}}    // namespace

