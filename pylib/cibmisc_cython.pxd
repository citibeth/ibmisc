from cibmisc cimport *

cdef extern from "ibmisc_cython.hpp" namespace "ibmisc::cython":
    cdef object linear_Weighted_shape(linear_Weighted &) except +
    cdef object linear_Weighted_type(linear_Weighted &) except +

    cdef object linear_Weighted_apply_weight(linear_Weighted &,
        int, PyObject *) except +

    cdef object linear_Weighted_apply_M(linear_Weighted &,
        PyObject *, double, bool) except +

    cdef object linear_Weighted_to_coo(linear_Weighted &) except +

    cdef object linear_Weighted_get_weights(linear_Weighted &, int) except +

    # Used for unit testing
    cdef unique_ptr[linear_Weighted] example_linear_weighted(string &) except +

