#
#  Source files and headers in this directory
#
file(GLOB HEOM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HEOM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HEOM_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(heom SHARED ${HEOM_SRC})
ADD_LIBRARY(heom_stat STATIC ${HEOM_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in heom")
   target_precompile_headers(heom REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(heom_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(heom      util_stat integrators_stat linalg_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(heom_stat util_stat integrators_stat linalg_stat specialfunctions_stat ${ext_libs})


