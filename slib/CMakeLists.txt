# ----------- ibmisc_c
set(IBMISC_SOURCE
    ibmisc/string.cpp
    ibmisc/datetime.cpp
    ibmisc/ibmisc.cpp
    ibmisc/geodesy.cpp
    ibmisc/VarTransformer.cpp
    ibmisc/ConstantSet.cpp
    ibmisc/ncfile.cpp
    ibmisc/indexing.cpp)


if (USE_PROJ4)
    list(APPEND IBMISC_SOURCE ibmisc/Proj2.cpp)
endif()

if (USE_BOOST)
    if (USE_NETCDF)
        list(APPEND IBMISC_SOURCE
            ibmisc/netcdf.cpp)
    endif()
endif()

if (USE_UDUNITS2)
    list(APPEND IBMISC_SOURCE
        ibmisc/udunits2.cpp)
endif()

if (USE_EVERYTRACE)
    add_library (ibmisc SHARED ${EVERYTRACE_c_REFADDR} ${IBMISC_SOURCE})
else()
    add_library(ibmisc SHARED ${IBMISC_SOURCE})
endif()

target_link_libraries(ibmisc ${EXTERNAL_LIBS})

install(TARGETS ibmisc DESTINATION lib)

# ------------ Header Files
## install Ibmisc headers
#file(GLOB IBMISC_HEADERS "ibmisc/*.hpp" )
#
## http://www.vtk.org/Wiki/CMake_FAQ#How_to_convert_a_semicolon_separated_list_to_a_whitespace_separated_string.3F
#foreach(arg ${IBMISC_HEADERS})
#  set(IBMISC_HEADERS_SPACES "${IBMISC_HEADERS_SPACES} ${arg}")
#endforeach(arg ${IBMISC_HEADERS})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")
install(FILES ${PROJECT_SOURCE_DIR}/slib/ibmisc/f90blitz_f.py DESTINATION include/ibmisc)


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


# Optional Items in libraries:
#USE_PROJ4
#  Proj2.cpp
#  Proj.hpp Proj2.hpp
#
#USE_BLITZ
#   f90blitz_f.f90
#  VarTransformer.hpp blitz.hpp f90blitz.hpp f90blitz_f.h
#
#USE_BOOST
#  boost/enum/*.hpp boost/enum.hpp
#
#USE_NETCDF && USE_BOOST
#  netcdf.cpp
#  netcdf.hpp
#
#USE_UDUNITS2
#  udunits2.hpp
#  udunits2.cpp
#
#
