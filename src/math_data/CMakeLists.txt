#
#  Source files and headers in this directory
#
file(GLOB DATA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB DATA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${DATA_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(data SHARED ${DATA_SRC})
ADD_LIBRARY(data_stat STATIC ${DATA_SRC})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_data")
   target_precompile_headers(data REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(data_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(data      ${ext_libs})
TARGET_LINK_LIBRARIES(data_stat ${ext_libs})



