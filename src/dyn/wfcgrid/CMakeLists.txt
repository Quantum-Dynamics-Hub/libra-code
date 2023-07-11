#
#  Source files and headers in this directory
#
file(GLOB WFCGRID_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB WFCGRID_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${WFCGRID_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(wfcgrid SHARED ${WFCGRID_SRC})
ADD_LIBRARY(wfcgrid_stat STATIC ${WFCGRID_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in wfcgrid")
   target_precompile_headers(wfcgrid REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(wfcgrid_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(wfcgrid      atomistic_stat linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(wfcgrid_stat atomistic_stat linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs})


