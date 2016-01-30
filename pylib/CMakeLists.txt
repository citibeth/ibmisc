
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

#		# http://stackoverflow.com/questions/12253389/numpy-c-api-link-several-object-files
#		set_property(
#			SOURCE ${CYTHON_MODULE}.cpp
#			PROPERTY COMPILE_DEFINITIONS USE_IMPORT
#		)


		# Remove NO_IMPORT_ARRAY compiler definition
		# http://stackoverflow.com/questions/24238937/how-to-change-a-compiler-flag-for-just-one-executable-in-cmake
#		get_source_file_property( MYPROPS ${CYTHON_MODULE}.cpp COMPILE_DEFINITIONS)
#message("xxxxxxxxxxxxxxxxxxxxxx " ${MYPROPS})
#		STRING( REPLACE "-DNO_IMPORT_ARRAY" "" MYPROPS ${MYPROPS} )
#		set_source_files_properties( ${CYTHON_MODULE}.cpp COMPILE_DEFINITIONS ${MYPROPS} )



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

	install(TARGETS ${CYTHON_MODULE}_so DESTINATION pylib)

endforeach()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/pylib/ DESTINATION pylib
          FILES_MATCHING PATTERN "*.py")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/pylib/ DESTINATION pylib
          FILES_MATCHING PATTERN "*.pxd")

include_directories(${PROJECT_SOURCE_DIR}/slib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

