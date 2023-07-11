#
#  Source files and headers in this directory
#
file(GLOB GRAPH_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB GRAPH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${GRAPH_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(graph      SHARED ${GRAPH_SRC})
ADD_LIBRARY(graph_stat STATIC ${GRAPH_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_data")
   target_precompile_headers(graph REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(graph_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(graph      ${ext_libs})
TARGET_LINK_LIBRARIES(graph_stat ${ext_libs})

                                                 
