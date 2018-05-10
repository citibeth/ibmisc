# IBMisc: Misc. Routines for IceBin (and other code)
# Copyright (c) 2013-2016 by Elizabeth Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
cimport numpy as np
np.import_array()
import scipy.sparse

cimport cibmisc            # Exports
cimport cibmisc_cython    # Internal
from cpython.object cimport *
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.memory cimport unique_ptr
from libcpp cimport bool
import functools
import operator

cibmisc.init()

cdef class NcIO:
    # cdef cibmisc.NcIO *cself    # Already defined in ibmisc.pxd

    def __cinit__(self, filePath, sfMode):
        self.cself = cibmisc.new_ncio(filePath.encode(), sfMode.encode())

    def __dealloc__(self):
        del self.cself

    def close(self):
        self.cself.close()

# blitz::Array --> Numpy Array
# http://stackoverflow.com/questions/23872946/force-numpy-ndarray-to-take-ownership-of-its-memory-in-cython

# Good Numpy / Cython tutorial
# https://github.com/cython/cython/wiki/tutorials-numpy


# This allows access to parts of PyArrayObject that the standard
# numpy.pxd does not.
# https://github.com/cython/cython/wiki/WrappingNumpy
#cdef extern from "numpy/arrayobject.h":
#    ctypedef int intp
#    ctypedef extern class numpy.ndarray [object PyArrayObject]:
#        cdef char *data
#        cdef int nd
#        cdef intp *dimensions
#        cdef intp *strides
#        cdef int flags

# ===================================================

cdef split_shape(ashape, alen):
    # See if we can find some set of dimensions matching ilen
    cumlen = np.cumprod(tuple(reversed(ashape)))
    try:
        icut = len(ashape) - 1 - next(i for i,x in enumerate(cumlen) if x==alen)
    except StopIteration:
        raise ValueError('Cannot find trailing dimension of {} in input of shape {}'.format(alen, ashape)) from None

    return ashape[0:icut], (
        functools.reduce(operator.mul, ashape[0:icut], 1),
        functools.reduce(operator.mul, ashape[icut:], 1))


cdef extern from "examples.hpp" namespace "ibmisc::cython":
    cdef void cyexample_double_blitz(PyObject *) except +
    cdef object cyexample_sparse_matrix() except +

# https://bitbucket.org/binet/cy-cxxfwk/src/c2dcf8bd9030b90fb59a168ebf293bb54ea4cb3f/cy/cyfwk.pyx?at=default&fileviewer=file-view-default
def example_double_blitz(a):
    cyexample_double_blitz(<PyObject *>a)

def example_sparse_matrix():
    data, shape = cyexample_sparse_matrix()
    return scipy.sparse.coo_matrix(data, shape)

# --------------------------
cdef class linear_Weighted:
    """The result of RegridMatrices.matrix()"""
    # cdef unique_ptr[cibmisc.linear_Weighted] cself    # Defined in ibmisc.pxd

    def __cinit__(self):
        self.cself = NULL

    def __dealloc__(self):
        if (self.cself):
            del self.cself

    @property
    def type(self):
        return cibmisc_cython.linear_Weighted_type(self.cself[0])

    @property
    def shape(self):
        return cibmisc_cython.linear_Weighted_shape(self.cself[0])

    def apply_weight(self, int dim, A_s):
        """Computes dot product of a weight vector with A_s.
        dim:
            0 = use weight vector for B (output) dimension
            1 = use weight vector for A (input) dimension
        A_s: Either:
            - A single vector (1-D array) to be transformed.
            - A 2-D array of row vectors to be transformed."""

        # Number of elements in sparse in put vector
        _,alen = self.shape
        leading_shape, new_shape = split_shape(A_s.shape, alen)
        A_s = A_s.reshape(new_shape)
        B_s = cibmisc_cython.linear_Weighted_apply_weight(self.cself[0], dim, <PyObject *>A_s)

        return B_s

    def apply_wM(self, A_s):
        return self.apply_weight(0,A_s)
    def apply_Mw(self, A_s):
        return self.apply_weight(1,A_s)



    def apply_M(self, A_s, fill=np.nan, bool force_conservation=True):
        """Applies the regrid matrix to A_s.
        A_s: Either:
            - A single vector (1-D array) to be transformed.
            - A 2-D array of row vectors to be transformed.
        fill:
            Un-set indices in output array will get this value.
        force_conservation: bool
            If M is not conservative, apply a conservation correction
            at the end"""


        # Number of elements in sparse in put vector
        _,alen = self.shape
        leading_shape, new_shape = split_shape(A_s.shape, alen)
        A_s = A_s.reshape(new_shape)
        B_s = cibmisc_cython.linear_Weighted_apply_M(self.cself[0],
            <PyObject *>A_s, fill, force_conservation)

        return B_s

    def to_coo(self):
        (data,shape) = cibmisc_cython.linear_Weighted_to_coo(self.cself[0])
        return scipy.sparse.coo_matrix(data, shape)

    def get_weights(self, int idim):
        return cibmisc_cython.linear_Weighted_get_weights(self.cself[0], idim)

    def ncio(self, NcIO ncio, vname):
        self.cself.ncio(ncio.cself[0], vname.encode())

# https://stackoverflow.com/questions/12204441/passing-c-pointer-as-argument-into-cython-function
cdef linear_Weighted_init(cibmisc.linear_Weighted *_cself):
    ret = linear_Weighted()
    ret.cself = _cself
    return ret

def example_linear_weighted(slinear_type):
    return linear_Weighted_init(
        cibmisc_cython.example_linear_weighted(slinear_type.encode()).release())

def nc_read_weighted(NcIO ncio, vname):
    return linear_Weighted_init(
        cibmisc.nc_read_weighted(ncio.cself[0].nc, vname.encode()).release())

