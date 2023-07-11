#
#  Source files and headers in this directory
#
file(GLOB OPERATORS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB OPERATORS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${OPERATORS_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(operators SHARED ${OPERATORS_SRC})
ADD_LIBRARY(operators_stat STATIC ${OPERATORS_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_operators")
   target_precompile_headers(operators REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(operators_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(operators      ${ext_libs})
TARGET_LINK_LIBRARIES(operators_stat ${ext_libs})


