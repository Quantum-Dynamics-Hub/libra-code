#
#  Source files and headers in this directory
#
file(GLOB IO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB IO_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${IO_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(io      SHARED ${IO_SRC})
ADD_LIBRARY(io_stat STATIC ${IO_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in io")
   target_precompile_headers(io REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(io_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(io      ${ext_libs})
TARGET_LINK_LIBRARIES(io_stat ${ext_libs})

                                                 
