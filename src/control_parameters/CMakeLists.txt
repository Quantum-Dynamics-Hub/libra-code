#
#  Source files and headers in this directory
#
file(GLOB CONTROL_PARAMETERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CONTROL_PARAMETERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CONTROL_PARAMETERS_HEADERS}) 



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(control_parameters SHARED ${CONTROL_PARAMETERS_SRC})
ADD_LIBRARY(control_parameters_stat STATIC ${CONTROL_PARAMETERS_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in io")
   target_precompile_headers(control_parameters REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(control_parameters_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( control_parameters      common_types linalg)
TARGET_LINK_LIBRARIES( control_parameters_stat common_types_stat linalg_stat)


