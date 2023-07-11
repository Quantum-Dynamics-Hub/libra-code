#
#  Source files and headers in this directory
#
file(GLOB GWP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB GWP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${GWP_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(gwp SHARED ${GWP_SRC})
ADD_LIBRARY(gwp_stat STATIC ${GWP_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in gwp")
   target_precompile_headers(gwp REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(gwp_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(gwp      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(gwp_stat linalg_stat ${ext_libs})


