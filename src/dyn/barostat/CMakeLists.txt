#
#  Source files and headers in this directory
#
file(GLOB BAROSTAT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB BAROSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${BAROSTAT_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(barostat SHARED ${BAROSTAT_SRC})
ADD_LIBRARY(barostat_stat STATIC ${BAROSTAT_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in barostat")
   target_precompile_headers(barostat REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(barostat_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(barostat      linalg_stat random_stat specialfunctions_stat io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(barostat_stat linalg_stat random_stat specialfunctions_stat io_stat ${ext_libs})


