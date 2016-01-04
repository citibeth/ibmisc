cmake_minimum_required(VERSION 3.1)
#enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if (NOT DEFINED USE_EVERYTRACE)
	set(USE_EVERYTRACE NO)
endif()

if (USE_EVERYTRACE)
	find_package(Everytrace REQUIRED)
	add_definitions(-DUSE_EVERYTRACE)
	include_directories(${EVERYTRACE_INCLUDE_DIR})
endif()

# -----------------------------------------------------
#https://cmake.org/pipermail/cmake/2007-February/012796.html
FIND_PACKAGE(Doxygen)

# ----------------------------------------------------------

include (setup_rpath)

# --------------------- Define our toolchain

find_package(gtest REQUIRED)
find_package(Blitz++ REQUIRED)
find_package(NetCDF4_C REQUIRED)
find_package(NetCDF4_CXX REQUIRED)

# -------- Process subdirectories of the build
add_subdirectory(slib)
#add_subdirectory(bin)
#add_subdirectory(test)

