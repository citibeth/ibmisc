# Input Variables
#    NETCDF_C_ROOT
# Produces:
#    NETCDF_C_LIBRARY
#    NETCDF_C_INCLUDE_DIR


FIND_PATH(NETCDF_C_INCLUDE_DIR netcdf.h
	HINTS ${NETCDF_C_ROOT}/include)

FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf
	HINTS ${NETCDF_C_ROOT}/lib)

IF (NETCDF_C_INCLUDE_DIR AND NETCDF_C_LIBRARY)
   SET(NETCDF_C_FOUND TRUE)
ENDIF (NETCDF_C_INCLUDE_DIR AND NETCDF_C_LIBRARY)

IF (NETCDF_C_FOUND)
   IF (NOT NETCDF_C_FIND_QUIETLY)
      MESSAGE(STATUS "Found NETCDF_C_LIBRARY: ${NETCDF_C_LIBRARY}")
   ENDIF (NOT NETCDF_C_FIND_QUIETLY)
ELSE (NETCDF_C_FOUND)
   IF (NETCDF_C_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find NETCDF_C")
   ENDIF (NETCDF_C_FIND_REQUIRED)
ENDIF (NETCDF_C_FOUND)
