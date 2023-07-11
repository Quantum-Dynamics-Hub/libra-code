#
#  Source files and headers in this directory
#
file(GLOB INTEGRATORS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB INTEGRATORS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${INTEGRATORS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(integrators      SHARED ${INTEGRATORS_SRC})
ADD_LIBRARY(integrators_stat STATIC ${INTEGRATORS_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: integrators")
   target_precompile_headers(integrators REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(integrators_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(integrators linalg meigen util specialfunctions)

TARGET_LINK_LIBRARIES(integrators_stat linalg_stat meigen_stat util_stat specialfunctions_stat)

