# Stuff defined in ibmisc.pyx
cimport cibmisc

cdef class NcIO:
	cdef cibmisc.NcIO *cself
