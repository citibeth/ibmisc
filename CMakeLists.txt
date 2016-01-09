cmake_minimum_required(VERSION 3.1)
#enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include_directories(${PROJECT_SOURCE_DIR}/slib)
#list (EXTERNAL_LIBS )

# -----------------------------------------
if (NOT DEFINED USE_EVERYTRACE)
	set(USE_EVERYTRACE NO)
endif()

if (USE_EVERYTRACE)
	find_package(Everytrace REQUIRED)
	add_definitions(-DUSE_EVERYTRACE)
	include_directories(${EVERYTRACE_INCLUDE_DIR})
	list(APPEND EXTERNAL_LIBS ${EVERYTRACE_LIBRARY})
endif()
# -----------------------------------------
if (NOT DEFINED USE_PROJ4)
	set(USE_PROJ4 YES)
endif()

if (USE_PROJ4)
	find_package(PROJ4 REQUIRED)
	add_definitions(-DUSE_PROJ4)
	include_directories(${PROJ4_INCLUDES})
	list(APPEND EXTERNAL_LIBS ${PROJ4_LIBRARIES})
endif()
# -----------------------------------------
if (NOT DEFINED USE_BLITZ)
	set(USE_BLITZ YES)
endif()
if (USE_BLITZ)
	find_package(blitz REQUIRED)
	add_definitions(-DUSE_BLITZ)
	include_directories(${BLITZ_INCLUDE_DIR})
	list(APPEND EXTERNAL_LIBS ${BLITZ_LIBRARY})
endif()
# -----------------------------------------
if (NOT DEFINED USE_NETCDF)
	set(USE_NETCDF YES)
endif()
if (USE_NETCDF)
	find_package(NetCDF4_C REQUIRED)
	find_package(NetCDF4_CXX REQUIRED)
	add_definitions(-DUSE_NETCDF)
	include_directories(
		${NETCDF4_C_INCLUDE_DIR}
		${NETCDF4_CXX_INCLUDE_DIR})
	list(APPEND EXTERNAL_LIBS ${NETCDF4_C_LIBRARY} ${NETCDF4_CXX_LIBRARY})
endif()
# -----------------------------------------------------
if (NOT DEFINED USE_BOOST)
	set(USE_BOOST YES)
endif()
if (USE_BOOST)
	find_package(boost COMPONENTS operators optional iterator preprocessor)
	add_definitions(-DUSE_BOOST)
	include_directories(${BOOST_INCLUDE_DIRS})
	list(APPEND EXTERNAL_LIBS ${BOOST_LIBRARIES})
endif()
# -----------------------------------------------------
if (NOT DEFINED USE_UDUNITS2)
	set(USE_UDUNITS2 YES)
endif()
if (USE_UDUNITS2)
	find_package(UDUNITS2 COMPONENTS operators optional iterator preprocessor)
	add_definitions(-DUSE_UDUNITS2)
	include_directories(${UDUNITS2_INCLUDES})
	list(APPEND EXTERNAL_LIBS ${UDUNITS2_LIBRARIES})
endif()
# -----------------------------------------------------
if (NOT DEFINED USE_GTEST)
	set(USE_GTEST YES)
endif()
# -----------------------------------------------------
#https://cmake.org/pipermail/cmake/2007-February/012796.html
FIND_PACKAGE(Doxygen)
# ----------------------------------------------------------

include (setup_rpath)


# -------- Process subdirectories of the build
add_subdirectory(slib)
#add_subdirectory(bin)
if (USE_GTEST)
	find_package(gtest REQUIRED)
	add_subdirectory(test)
endif()
