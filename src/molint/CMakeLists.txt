#
#  Source files and headers in this directory
#
file(GLOB MOLINT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MOLINT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MOLINT_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(molint SHARED ${MOLINT_SRC})
ADD_LIBRARY(molint_stat STATIC ${MOLINT_SRC})



#if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
#   message(STATUS "Compiling using pre-compiled header support: molint")
#   target_precompile_headers(molint REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
#   target_precompile_headers(molint_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
#endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(molint      linalg specialfunctions)
TARGET_LINK_LIBRARIES(molint_stat linalg_stat specialfunctions_stat)


