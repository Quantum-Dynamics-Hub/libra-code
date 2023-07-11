#
#  Source files and headers in this directory
#
file(GLOB RANDOM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB RANDOM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${RANDOM_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(random SHARED ${RANDOM_SRC})
ADD_LIBRARY(random_stat STATIC ${RANDOM_SRC})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in random")
   target_precompile_headers(random REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(random_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(random      ${ext_libs})
TARGET_LINK_LIBRARIES(random_stat ${ext_libs})


