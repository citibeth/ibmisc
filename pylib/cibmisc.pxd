from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cblitz
from cpython.object cimport *

# This will make the C++ class def for Rectangle available..
cdef extern from "ibmisc/indexing.hpp" namespace "ibmisc":
	pass

	cdef cppclass Indexing[TupleT, IndexT]:
		vector[TupleT] base
		vector[TupleT] extent
		vector[int] indices

cdef extern from "ibmisc/netcdf.hpp" namespace "ibmisc":
	pass

	cdef cppclass NcIO:
		NcIO(string fname, int fMode) except +
		void close()

cdef extern from "ibmisc/cython.hpp" namespace "ibmisc::cython":
	cdef extern void init()
	cdef extern call_import_array()
	cdef NcIO *new_ncio(string fname, string sFileMode)
