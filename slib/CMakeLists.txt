include_directories(${PROJECT_SOURCE_DIR}/slib)


# ----------- ibmisc_c
set(IBMISC_SOURCE ibmisc/ncutil.cpp)

add_library(ibmisc SHARED
	${EVERYTRACE_REFADDR}	# Blank if not USE_EVERYTRACE
	${IBMISC_SOURCE})

target_link_libraries(ibmisc ${EVERYTRACE_C_LIBRARY})
install(TARGETS ibmisc DESTINATION lib)

# ------------ Header Files
# install Ibmisc headers
file(GLOB IBMISC_HEADERS "giss/*.hpp")

# http://www.vtk.org/Wiki/CMake_FAQ#How_to_convert_a_semicolon_separated_list_to_a_whitespace_separated_string.3F
foreach(arg ${IBMISC_HEADERS})
  set(IBMISC_HEADERS_SPACES "${IBMISC_HEADERS_SPACES} ${arg}")
endforeach(arg ${IBMISC_HEADERS})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")


# ------------- API Docs
IF (DOXYGEN_FOUND)

	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

	add_custom_target(apidoc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)


#ADD_CUSTOM_TARGET(apidoc_forced
#  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
#  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
#  COMMAND ${CMAKE_COMMAND} -E echo "Done."
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ENDIF (DOXYGEN_FOUND)
