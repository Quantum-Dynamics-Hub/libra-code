#
#  Source files and headers in this directory
#
file(GLOB CELL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CELL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CELL_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(cell SHARED ${CELL_SRC})
ADD_LIBRARY(cell_stat STATIC ${CELL_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: cell")
   target_precompile_headers(cell REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(cell_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(cell      linalg)
TARGET_LINK_LIBRARIES(cell_stat linalg_stat)


