#
#  Source files and headers in this directory
#
file(GLOB OPT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB OPT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${OPT_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(opt      SHARED ${OPT_SRC})
ADD_LIBRARY(opt_stat STATIC ${OPT_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: opt")
   target_precompile_headers(opt REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(opt_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(opt       linalg  io ${ext_libs} )
TARGET_LINK_LIBRARIES(opt_stat  linalg_stat  io_stat ${ext_libs} )




