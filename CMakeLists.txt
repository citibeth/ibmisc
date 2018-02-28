cmake_minimum_required(VERSION 3.1)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include_directories(${PROJECT_SOURCE_DIR}/slib)
#list (EXTERNAL_LIBS )

# --- Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
list(APPEND EXTERNAL_LIBS ${ZLIB_LIBRARIES})

# -----------------------------------------
find_package(Everytrace REQUIRED)
include_directories(${EVERYTRACE_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBS ${EVERYTRACE_LIBRARY})
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
    find_package(Blitz REQUIRED)
    add_definitions(-DUSE_BLITZ)
    include_directories(${BLITZ_INCLUDE_DIR})
    list(APPEND EXTERNAL_LIBS ${BLITZ_LIBRARY})
endif()
# -----------------------------------------
if (NOT DEFINED USE_NETCDF)
    set(USE_NETCDF YES)
endif()
if (USE_NETCDF)
#   find_package(NetCDF4_C REQUIRED)
    find_package(NetCDF4_CXX REQUIRED)
    add_definitions(-DUSE_NETCDF)
    include_directories(
#       ${NETCDF4_C_INCLUDE_DIR}
        ${NETCDF4_CXX_INCLUDE_DIR})
#   list(APPEND EXTERNAL_LIBS ${NETCDF4_C_LIBRARY} ${NETCDF4_CXX_LIBRARY})
    list(APPEND EXTERNAL_LIBS ${NETCDF4_CXX_LIBRARY})
endif()
# -----------------------------------------------------
if (NOT DEFINED USE_BOOST)
    set(USE_BOOST YES)
endif()
if (USE_BOOST)
    # We USE these parts of Boost... but they aren't components in the CMake
    # sense because they are header-only:
    # find_package(Boost COMPONENTS operators optional iterator preprocessor filesystem)
    find_package(Boost COMPONENTS filesystem)
    add_definitions(-DUSE_BOOST)
    include_directories(${Boost_INCLUDE_DIRS})
    list(APPEND EXTERNAL_LIBS ${Boost_LIBRARIES})
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
if (NOT DEFINED BUILD_PYTHON)
    set(BUILD_PYTHON YES)
endif()
# -----------------------------------------------------
if (NOT DEFINED USE_GTEST)
    set(USE_GTEST YES)
endif()
# -----------------------------------------------------
#https://cmake.org/pipermail/cmake/2007-February/012796.html
if (NOT DEFINED BUILD_DOCS)
    set(BUILD_DOCS NO)
endif()
if (BUILD_DOCS)
    FIND_PACKAGE(Doxygen)
endif()
# ----------------------------------------------------------

include (setup_rpath)

# Include all the transitive dependencies determined by EasyBuild.
# If we're not running with EasyBuild, this does nothing...
include_directories($ENV{SPACK_TRANSITIVE_INCLUDE_PATH})


# -------- Process subdirectories of the build
add_subdirectory(slib)
#add_subdirectory(bin)
if (USE_GTEST)
    find_package(Gtest REQUIRED)
    add_subdirectory(tests)
endif()

if (BUILD_PYTHON)
    add_subdirectory(pylib)
endif()
