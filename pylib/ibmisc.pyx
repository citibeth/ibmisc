import numpy as np
cimport numpy as np
np.import_array()
import scipy.sparse

cimport cibmisc
from cpython.object cimport *
from cython.operator cimport dereference as deref, preincrement as inc

cibmisc.init()

cdef class NcIO:
	# cdef cibmisc.NcIO *cself

	def __cinit__(self, filePath, sfMode):
		self.cself = cibmisc.new_ncio(filePath.encode(), sfMode.encode());

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
