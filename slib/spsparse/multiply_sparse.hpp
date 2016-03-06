#ifndef SPSPARSE_MULTIPLY_SPARSE_HPP
#define SPSPARSE_MULTIPLY_SPARSE_HPP

#include <spsparse/xiter.hpp>
#include <spsparse/array.hpp>

namespace spsparse {

/** @defgroup multiply_sparse multiply_sparse.hpp
@brief Fully sparse-sparse multiplication algorithms.

Code Example
@code
// --------- Matrix-Matrix Multiplation
VectorCooMatrix<int, double> A,B,AB,BtAt;
multiply(AB, 1.0, NULL, A,'.', NULL, B,'.', NULL);
multiply(BtAt, 1.0, NULL, B,'T', NULL, A,'T', NULL);
// Now: AB == BtAt (but we don't yet have a way to compare arrays)

VectorCooMatrix<int, double> ret;
VectorCooVector<int, double> scalei(A.shape[0]);
VectorCooVector<int, double> scalej(A.shape[1]);
multiply(AB, 17.0, &scale, A,'.', &scalej, B,'.', NULL);

// --------- Matrix-Vector Multiplation
VectorCooVector<int, double> V(B.shape[1]);
VectorCooVector<int, double> retV;
multiply(retV, 1.0, NULL, AB,'.', NULL, V);
@endcode


@{
*/

// -------------------------------------------------------------

/** @brief Internal helper class for sparse-sparse multiplication.

Superclass used to chose whether we do/do not multiply our
matrix by a diagonal scaling matrix */
template<class MatT>
struct MultXiter {
    virtual typename MatT::index_type index() = 0;
    virtual bool eof() = 0;
    virtual void operator++() = 0;
    virtual typename MatT::val_type scale_val() = 0;
    virtual typename DimBeginningsXiter<MatT>::sub_xiter_type sub_xiter() = 0;
};


/** @brief Internal helper class for sparse-sparse multiplication.

Use this if you do NOT want to use a diagonal scaling matrix. */
template<class MatT>
class SimpleMultXiter : public MultXiter<MatT>
{
    // SPSPARSE_LOCAL_TYPES(MatT);
    DimBeginningsXiter<MatT> ii;
public:
    SimpleMultXiter(MatT const &A)
        : ii(A.dim_beginnings_xiter())
    {}

    typename MatT::index_type index() { return *ii; }
    bool eof() { return ii.eof(); }
    void operator++() { ++ii; }
    typename MatT::val_type scale_val() { return 1; }
    typename DimBeginningsXiter<MatT>::sub_xiter_type sub_xiter() { return ii.sub_xiter(); }
};
// ---------------------------------------------------------
/** @brief Internal helper class for sparse-sparse multiplication.

Use this if you DO want to use a diagonal scaling matrix. */
template<class MatT, class ScaleT>
class ScaledMultXiter : public MultXiter<MatT>
{
    // SPSPARSE_LOCAL_TYPES(MatT);

    Join2Xiter<
        DimBeginningsXiter<MatT>,
        ValSTLXiter<typename ScaleT::const_dim_iterator>> ii;
public:
    ScaledMultXiter(MatT const &A, ScaleT const &scale)
        : ii(A.dim_beginnings_xiter(),
            make_val_xiter(scale.dim_begin(0), scale.dim_end(0)))
    {}

    typename MatT::index_type index() { return *ii.i1; }
    bool eof() { return ii.eof(); }
    void operator++() { ++ii; }
    typename MatT::val_type scale_val() { return ii.i2.val(); }
    typename DimBeginningsXiter<MatT>::sub_xiter_type sub_xiter() { return ii.i1.sub_xiter(); }
};

// ---------------------------------------------------------
/** @brief Internal helper function for sparse-sparse multiplication.

Chooses the correct spsparse::MultXiter, based on whether or not there
is a diagonal scaling matrix.
*/
template<class MatT, class ScaleT>
std::unique_ptr<MultXiter<MatT>> new_mult_xiter(MatT const &A, ScaleT const *scale);

template<class MatT, class ScaleT>
std::unique_ptr<MultXiter<MatT>> new_mult_xiter(MatT const &A, ScaleT const *scale)
{
    if (scale) {
        return std::unique_ptr<MultXiter<MatT>>(
            new ScaledMultXiter<MatT, ScaleT>(A, *scale));
    } else {
        return std::unique_ptr<MultXiter<MatT>>(
            new SimpleMultXiter<MatT>(A));
    }

}
// -------------------------------------------------------------
/** @brief Matrix-matrix multiply.

Computes:
@code
ret = C * diag(scalei) * A * diag(scalej) * B * diag(scalek)
@endcode

@param ret Accumulator for output.
@param C Constant to multiply by.
@param scalei Vector for diagonal scale matrix (NULL if not used).
@param A Matrix
@param transpose_A If =='T', then use transpose(A) instead of A
@param scalej Vector for diagonal scale matrix (NULL if not used).
@param B Matrix
@param transpose_B If =='T', then use transpose(B) instead of A
@param scalek Vector for diagonal scale matrix (NULL if not used).
@param duplicate_policy Use if A or B needs to be consolidated.
@param zero_nan Use if A or B needs to be consolidated.

@see spsparse::consolidate()
*/
template<class ScaleIT, class MatAT, class ScaleJT, class MatBT, class ScaleKT, class AccumulatorT>
void multiply(
    AccumulatorT &ret,
    double C,       // Multiply everything by this
    ScaleIT const *scalei,
    MatAT const &A,
    char transpose_A,           // 'T' for transpose, '.' otherwise
    ScaleJT const *scalej,      // Vector
    MatBT const &B,
    char transpose_B,           // 'T' for transpose, '.' otherwise
    ScaleKT const *scalek,
    DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
    bool zero_nan = false);

template<class ScaleIT, class MatAT, class ScaleJT, class MatBT, class ScaleKT, class AccumulatorT>
void multiply(
    AccumulatorT &ret,
    double C,       // Multiply everything by this
    ScaleIT const *scalei,
    MatAT const &A,
    char transpose_A,           // 'T' for transpose, '.' otherwise
    ScaleJT const *scalej,      // Vector
    MatBT const &B,
    char transpose_B,           // 'T' for transpose, '.' otherwise
    ScaleKT const *scalek,
    DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
    bool zero_nan = false)
{
    std::array<int,2> const &adims = (transpose_A == 'T' ? COL_MAJOR : ROW_MAJOR);
    std::array<int,2> const &bdims = (transpose_B == 'T' ? COL_MAJOR : ROW_MAJOR);


    // Set dimensions of output, even if we store nothing in it.
    std::array<int,2> const &a_sort_order(transpose_A == 'T' ? COL_MAJOR : ROW_MAJOR);
    std::array<int,2> const &b_sort_order(transpose_B == 'T' ? ROW_MAJOR : COL_MAJOR);
    ret.set_shape({A.shape[adims[0]], B.shape[bdims[1]]});

    // Check inner dimensions
    if (A.shape[adims[1]] != B.shape[bdims[0]]) {
        (*spsparse_error)(-1, "Inner dimensions for A (%ld) and B (%ld) must match!", A.shape[adims[1]], B.shape[bdims[0]]);
    }

    if (scalei && A.shape[adims[0]] != scalei->shape[0]) {
        (*spsparse_error)(-1, "Dimension for scalei (%ld) must match A.shape[0]=%ld\n", scalei->shape[0], A.shape[adims[0]]);
    }

    if (scalej && A.shape[adims[1]] != scalej->shape[0]) {
        (*spsparse_error)(-1, "Dimension for scalej (%ld) must match A.shape[1]=%ld\n", scalej->shape[0], A.shape[adims[1]]);
    }
    if (scalek && B.shape[bdims[1]] != scalek->shape[0]) {
        (*spsparse_error)(-1, "Dimension for scalek (%ld) must match B.shape[1]=%ld\n", scalek->shape[0], B.shape[bdims[1]]);
    }


    // Short-circuit return on empty output
    // (DimBeginningXiter doesn't like size()==0)
    if ((isnone(C)) 
        || (scalei && scalei->size() == 0)
        || (A.size() == 0)
        || (scalej && scalej->size() == 0)
        || (B.size() == 0)
        || (scalek && scalek->size() == 0))
    { return; }

    // --------- Consolidate the matrices if needed
    Consolidate<MatAT> Acon(&A, a_sort_order, duplicate_policy, zero_nan);
    Consolidate<MatBT> Bcon(&B, b_sort_order, duplicate_policy, zero_nan);

    // Multiply each row by each column
    // ------ Loop 1: Rows in A 
    for (auto join_a(new_mult_xiter(Acon(), scalei));
        !join_a->eof(); ++(*join_a))
    {
        if (isnone(join_a->scale_val())) continue;
        auto aix(join_a->index());
        auto a_scale(join_a->scale_val());

#if 0
printf("Starting row %d:", aix);
for (auto ii=join_a->sub_xiter(); !ii.eof(); ++ii) {
    printf(" (%d : %g)", *ii, ii.val());
    }
printf("\n");
#endif

        // ------ Loop 2: Cols in B
        for (auto join_b(new_mult_xiter(Bcon(), scalek));
            !join_b->eof(); ++(*join_b))
        {
            if (isnone(join_b->scale_val())) continue;
            auto bix(join_b->index());
            auto b_scale(join_b->scale_val());
#if 0
printf("  Starting col %d\n", bix);
#endif

            // ---------- Loop 3: Multiply A row by B column
            typename AccumulatorT::val_type sum = 0;

            if (scalej) {
                // ----- With scalej
                for (auto ab(join3_xiter(
                    join_a->sub_xiter(),
                    make_val_xiter(scalej->dim_begin(0), scalej->dim_end(0)),
                    join_b->sub_xiter()));
                    !ab.eof(); ++ab)
                { sum += ab.i1.val() * ab.i2.val() * ab.i3.val(); }
            } else {
                // ----- Without scalej
                for (auto ab(join2_xiter(
                    join_a->sub_xiter(),
                    join_b->sub_xiter()));
                    !ab.eof(); ++ab)
                { sum += ab.i1.val() * ab.i2.val(); }
            }

            if (!isnone(sum)) {
#if 0
printf("set: (%d %d : %g)\n", aix, bix, sum * C * a_scale * b_scale);
#endif
                ret.add({aix, bix}, sum * C * a_scale * b_scale);
            }

        }
    }

}
// ------------------------------------------------------------
/** @brief Matrix-vector multiply.

Computes:
@code
ret = C * diag(scalei) * A * diag(scalej) * V
@endcode

@param ret Accumulator for output.
@param C Constant to multiply by.
@param scalei Vector for diagonal scale matrix (NULL if not used).
@param A Matrix
@param transpose_A If =='T', then use transpose(A) instead of A
@param scalej Vector for diagonal scale matrix (NULL if not used).
@param V Vector
@param duplicate_policy Use if A or B needs to be consolidated.
@param zero_nan Use if A or B needs to be consolidated.

@see spsparse::consolidate()
*/
template<class ScaleIT, class MatAT, class ScaleJT, class VecT, class AccumulatorT>
void multiply(
    AccumulatorT &ret,
    double C,       // Multiply everything by this
    ScaleIT const *scalei,
    MatAT const &A,
    char transpose_A,           // 'T' for transpose, '.' otherwise
    ScaleJT const *scalej,      // Vector
    VecT const &V,
    DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
    bool zero_nan = false);

template<class ScaleIT, class MatAT, class ScaleJT, class VecT, class AccumulatorT>
void multiply(
    AccumulatorT &ret,
    double C,       // Multiply everything by this
    ScaleIT const *scalei,
    MatAT const &A,
    char transpose_A,           // 'T' for transpose, '.' otherwise
    ScaleJT const *scalej,      // Vector
    VecT const &V,
    DuplicatePolicy duplicate_policy = DuplicatePolicy::ADD,
    bool zero_nan = false)
{
    // Set dimensions of output, even if we store nothing in it.
    std::array<int,2> const &a_sort_order(transpose_A == 'T' ? COL_MAJOR : ROW_MAJOR);
    ret.set_shape({A.shape[a_sort_order[0]]});

    // Check inner dimensions
    if (A.shape[a_sort_order[1]] != V.shape[0]) {
        (*spsparse_error)(-1, "Inner dimensions for A (%ld) and V (%ld) must match!", A.shape[a_sort_order[1]], V.shape[0]);
    }

    // Short-circuit return on empty output
    // (DimBeginningXiter doesn't like size()==0)
    if ((isnone(C)) 
        || (scalei && scalei->size() == 0)
        || (A.size() == 0)
        || (scalej && scalej->size() == 0)
        || (V.size() == 0))
    { return; }

    // --------- Consolidate the matrices if needed
    Consolidate<MatAT> Acon(&A, a_sort_order, duplicate_policy, zero_nan);
    Consolidate<VecT> Vcon(&V, {0}, duplicate_policy, zero_nan);

//std::cout << "Vcon: " << Vcon() << std::endl;

    // Multiply each row by each column
    // ------ Loop 1: Rows in A 
    for (auto join_a(new_mult_xiter(Acon(), scalei));
        !join_a->eof(); ++(*join_a))
    {
        if (isnone(join_a->scale_val())) continue;
        auto aix(join_a->index());
        auto a_scale(join_a->scale_val());

#if 0
printf("Starting row %d:", aix);
for (auto ii=join_a->sub_xiter(); !ii.eof(); ++ii) {
    printf(" (%d : %g)", *ii, ii.val());
    }
printf("\n");
#endif
        // ---------- Loop 3: Multiply A row by the single B column
        typename AccumulatorT::val_type sum = 0;

        if (scalej) {
            // ----- With scalej
            for (auto ab(join3_xiter(
                join_a->sub_xiter(),
                make_val_xiter(scalej->dim_begin(0), scalej->dim_end(0)),
                make_val_xiter(Vcon().dim_begin(0), Vcon().dim_end(0))));
                !ab.eof(); ++ab)
            { sum += ab.i1.val() * ab.i2.val() * ab.i3.val(); }
        } else {
            // ----- Without scalej
            for (auto ab(join2_xiter(
                join_a->sub_xiter(),
                make_val_xiter(Vcon().dim_begin(0), Vcon().dim_end(0))));
                !ab.eof(); ++ab)
            {
//printf("    triplet %d=?%d: %g*%g = %g\n", *ab.i1, *ab.i2, ab.i1.val(), ab.i2.val(), ab.i1.val()*ab.i2.val());
                sum += ab.i1.val() * ab.i2.val();
            }
        }

        if (!isnone(sum)) {
#if 0
printf("    set: (%d : %g)\n", aix, sum * C * a_scale);
#endif
            ret.add({aix}, sum * C * a_scale);
        }

    }

}

// ------------------------------------------------------------
template<class ScaleAT, class ScaleBT, class AccumulatorT>
void multiply_ele(
    AccumulatorT &ret,
    ScaleAT const &A,
    ScaleBT const &B);

template<class ScaleAT, class ScaleBT, class AccumulatorT>
void multiply_ele(
    AccumulatorT &ret,
    ScaleAT const &A,
    ScaleBT const &B)
{
    for (auto jii = join2_xiter(
        make_val_xiter(A.begin(), A.end()),
        make_val_xiter(B.begin(), B.end()));
        !jii.eof(); ++jii)
    {
        ret.add(*jii.i1, jii.i1.val() * jii.i2.val());
    }
}

/** @} */

} // namespace


#endif  // Guard
