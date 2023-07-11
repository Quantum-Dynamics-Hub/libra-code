#
#  Source files and headers in this directory
#
file(GLOB SYMMETRY_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SYMMETRY_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SYMMETRY_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(symmetry SHARED ${SYMMETRY_SRC})
ADD_LIBRARY(symmetry_stat STATIC ${SYMMETRY_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_symmetry")
   target_precompile_headers(symmetry REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(symmetry_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(symmetry      linalg ${ext_libs})
TARGET_LINK_LIBRARIES(symmetry_stat linalg_stat ${ext_libs})


