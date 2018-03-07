/** Return value of a sparse matrix */
struct WeightedSparse {
    ibmisc::TmpAlloc tmp;            // Stuff that needs same lifetime as WeightedSparse

    std::array<SparseSetT *,2> dims;            // Dense-to-sparse mapping for the dimensions

    // If M=BvA, then wM = wBvA = area of B cells
    DenseArrayT<1> wM;           // Dense indexing

    std::unique_ptr<EigenSparseMatrixT> M;    // Dense indexing

    // Area of A cells
    DenseArrayT<1> Mw;


    /** True if this regridding matrix is conservative.  Matrices could be
    non-conservative, for example, in the face of smoothing on I.  Or when
    regridding between the IceBin and ModelE ice sheets. */
    bool conservative;

    WeightedSparse(std::array<SparseSetT *,2> _dims, bool _conservative) : dims(_dims), conservative(_conservative) {}

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
    EigenDenseMatrixT apply_e(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable (_b means dense here)
        double fill = std::numeric_limits<double>::quiet_NaN(),    // Fill value for cells not in BvA matrix
        bool force_conservation=true) const;     // Set if you want apply_e() to conserve, even if !M->conservative


    /** Apply to multiple variables
    @return Blitz type */
    blitz::Array<double,2> apply(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
        double fill,    // Fill value for cells not in BvA matrix
        bool force_conservation,
        ibmisc::TmpAlloc &tmp) const;


    /** Apply to a single variable */
    blitz::Array<double,1> apply(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,1> const &A_b,       // A_b{j} One variable
        double fill,    // Fill value for cells not in BvA matrix
        bool force_conservation,
        ibmisc::TmpAlloc &tmp) const;

    /** Read/write to NetCDF */
    void ncio(ibmisc::NcIO &ncio,
        std::string const &vname,
        std::array<std::string,2> dim_names);

};
