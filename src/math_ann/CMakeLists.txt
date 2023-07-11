#
#  Source files and headers in this directory
#
file(GLOB ANN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ANN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${ANN_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(ann SHARED      ${ANN_SRC})
ADD_LIBRARY(ann_stat STATIC ${ANN_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_ann")
   target_precompile_headers(ann REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(ann_stat REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#

TARGET_LINK_LIBRARIES(ann      specialfunctions data linalg ${ext_libs} )
TARGET_LINK_LIBRARIES(ann_stat specialfunctions_stat data_stat linalg_stat ${ext_libs} )


