#
#  Source files and headers in this directory
#
file(GLOB MC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MC_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(montecarlo      SHARED ${MC_SRC})
ADD_LIBRARY(montecarlo_stat STATIC ${MC_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: montecarlo")
   target_precompile_headers(montecarlo REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(montecarlo_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(montecarlo      random  linalg  io ${ext_libs} )
TARGET_LINK_LIBRARIES(montecarlo_stat random_stat  linalg_stat  io_stat ${ext_libs} )




