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

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, unique_ptr
cimport cblitz
from cpython.object cimport *

#cdef extern from "<utility>" namespace "std":
#    cdef T std_move "std::move" [T](T t)

## See: https://github.com/apache/arrow/commit/000e1e34d6ad3b6c1a1bc430974f2eac05f96173
#cdef extern from "<memory>" namespace "std" nogil:
#
#    cdef cppclass unique_ptr[T]:
#        unique_ptr()
#        unique_ptr(T*)
#        T* get()
#        T* release()
#        void reset()
#        void reset(nullptr_t)
#        void reset(T*)
#        void swap(unique_ptr&)
#
#    cdef cppclass shared_ptr[T]:
#        shared_ptr()
#        shared_ptr(T*)
#        T* get()
#        void reset()
#        void reset(T* p)
#        void swap(shared_ptr&)
#
#cdef extern from "<array>" namespace "std" nogil:
#    cdef cppclass array[T,rank]:
#        array() except+
#        T &operator[](size_t)
#
# ==========================================================

# This will make the C++ class def for Rectangle available..
cdef extern from "ibmisc/indexing.hpp" namespace "ibmisc":
    pass

    cdef cppclass Indexing[TupleT, IndexT]:
        vector[TupleT] base
        vector[TupleT] extent
        vector[int] indices

cdef extern from "<netcdf>" namespace "netCDF":
    cdef cppclass NcGroup:
        pass

cdef extern from "ibmisc/netcdf.hpp" namespace "ibmisc":
    cdef cppclass NcIO:
        NcGroup *nc
        NcIO(string fname, int fMode) except +
        void close()


cdef extern from "ibmisc/cython.hpp" namespace "ibmisc::cython":
    cdef extern void init()
    cdef extern call_import_array()
    cdef NcIO *new_ncio(string fname, string sFileMode)


cdef extern from "<ibmisc/linear/linear.hpp>" namespace "ibmisc::linear":
    cdef cppclass linear_Weighted "ibmisc::linear::Weighted":
        void ncio(NcIO &, string &) except +

    cdef extern unique_ptr[linear_Weighted] nc_read_weighted(
        NcGroup *nc, string vname) except +


