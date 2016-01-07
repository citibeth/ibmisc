# https://cmake.org/cmake/help/v3.0/module/FindGTest.html

# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/slib)

SET(ALL_LIBS ${GTEST_BOTH_LIBRARIES} ${EVERYTRACE_C_LIBRARY} ibmisc)

add_executable(test_netcdf test_netcdf.cpp)
target_link_libraries(test_netcdf ${ALL_LIBS})
add_test(AllTests test_netcdf)

#add_executable(test_array test_array.cpp)
#target_link_libraries(test_array ${ALL_LIBS})
#add_test(AllTests test_array)
#
#add_executable(test_multiply_sparse test_multiply_sparse.cpp)
#target_link_libraries(test_multiply_sparse ${ALL_LIBS})
#add_test(AllTests test_multiply_sparse)

