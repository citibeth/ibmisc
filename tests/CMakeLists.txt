# ---------------------------------
# https://cmake.org/cmake/help/v3.0/module/FindGTest.html
# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/slib)


SET(ALL_LIBS ${GTEST_LIBRARY} ${EXTERNAL_LIBS} ibmisc)


foreach(TEST netcdf iter blitz indexing memory)
	add_executable(ibmisc_${TEST} ibmisc/test_${TEST}.cpp)
	target_link_libraries(ibmisc_${TEST} ${ALL_LIBS})
	add_test(AllTests ibmisc_${TEST})
endforeach()

foreach(TEST xiter array netcdf multiply_sparse)
	add_executable(spsparse_${TEST} spsparse/test_${TEST}.cpp)
	target_link_libraries(spsparse_${TEST} ${ALL_LIBS})
	add_test(AllTests spsparse_${TEST})
endforeach()
