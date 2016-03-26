# IBMisc: Misc. Routines for IceBin (and other code)
# Copyright (c) 2013-2016 by Elizabeth Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
