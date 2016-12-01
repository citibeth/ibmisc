find_package(Python REQUIRED)
find_package(Cython REQUIRED)
find_package(Numpy REQUIRED)

# ----------Stuff to help with running Python before this project is installed.

# --- Set up PYTHONPATH to find stuff
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/env.in
    ${CMAKE_CURRENT_BINARY_DIR}/env @ONLY)
execute_process(COMMAND chmod a+x ${CMAKE_CURRENT_BINARY_DIR}/env)

## --- Convenience link to Python source
#execute_process(COMMAND ln -s -h -f ${CMAKE_CURRENT_SOURCE_DIR} src
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# -------------------------------------------------------------

include_directories(${PYTHON_INCLUDES} ${NUMPY_INCLUDE_DIRS})
list(APPEND EXTERNAL_LIBS ${PYTHON_LIBRARY})
# Compile everything (including Cython extension) as a single
# "compilation unit" for Numpy.  See:
# http://stackoverflow.com/questions/12253389/numpy-c-api-link-several-object-files
add_definitions(-DPY_ARRAY_UNIQUE_SYMBOL=ibmisc_array_unique_symbol)

# Only ONE compilation unit should NOT have -DNO_IMPORT_ARRAY
foreach (SOURCE_FILE cython.cpp examples.cpp)
    list(APPEND PYLIB_SOURCES ${SOURCE_FILE})
    set_property(SOURCE ${SOURCE_FILE}
        PROPERTY COMPILE_DEFINITIONS NO_IMPORT_ARRAY )
endforeach()

foreach (CYTHON_MODULE ibmisc)
    add_custom_command(OUTPUT ${CYTHON_MODULE}.cpp
        MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/pylib/${CYTHON_MODULE}.pyx
        DEPENDS cibmisc.pxd
        COMMAND ${CYTHON_EXECUTABLE} --cplus -3
            --include-dir ${PROJECT_SOURCE_DIR}/pylib
            ${PROJECT_SOURCE_DIR}/pylib/${CYTHON_MODULE}.pyx
            -o ${CYTHON_MODULE}.cpp)

    ADD_LIBRARY(${CYTHON_MODULE}_so MODULE ${CYTHON_MODULE}.cpp ${PYLIB_SOURCES} ${EVERYTRACE_REFADDR})

    target_link_libraries(${CYTHON_MODULE}_so ibmisc ${PYTHON_LIBRARY})

    set_target_properties(${CYTHON_MODULE}_so
        PROPERTIES
        OUTPUT_NAME ${CYTHON_MODULE}
        PREFIX ""
        SUFFIX ".so"
        OUTPUT ${CYTHON_MODULE}.so)



    # Set RPATH in the installed library
    # http://www.cmake.org/Wiki/CMake_RPATH_handling
    # http://www.kitware.com/blog/home/post/510
    set_target_properties(${CYTHON_MODULE}_so PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE)

    install(TARGETS ${CYTHON_MODULE}_so DESTINATION lib/python${PYTHON_VERSION}/site-packages)

endforeach()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/pylib/ DESTINATION lib/python${PYTHON_VERSION}/site-packages
          FILES_MATCHING PATTERN "*.py")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/pylib/ DESTINATION include
          FILES_MATCHING PATTERN "*.pxd")

include_directories(${PROJECT_SOURCE_DIR}/slib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

