# ---------------------------------
# https://cmake.org/cmake/help/v3.0/module/FindGTest.html
# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(Gtest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/slib)


SET(ALL_LIBS ${GTEST_LIBRARY} ${EXTERNAL_LIBS} ibmisc)


foreach(TEST netcdf iter blitz indexing memory var_transformer constant_set datetime string filesystem bundle permutation runlength)
    add_executable(ibmisc_${TEST} ibmisc/test_${TEST}.cpp)
    target_link_libraries(ibmisc_${TEST} ${ALL_LIBS})
    add_test(AllTests ibmisc_${TEST})
endforeach()

foreach(TEST array netcdf)
    add_executable(spsparse_${TEST} spsparse/test_${TEST}.cpp)
    target_link_libraries(spsparse_${TEST} ${ALL_LIBS})
    add_test(AllTests spsparse_${TEST})
endforeach()

# This test has a second Fortran file in it
add_executable(ibmisc_fortranio ibmisc/test_fortranio.cpp ibmisc/help_fortranio.F90)
#add_executable(ibmisc_fortranio ibmisc/test_fortranio.cpp)
target_link_libraries(ibmisc_fortranio ${ALL_LIBS} gfortran)   # Linking in gfortran is GCC-only
add_test(AllTests ibmisc_fortranio)
