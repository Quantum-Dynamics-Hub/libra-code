#
#  Source files and headers in this directory
#
file(GLOB POT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB POT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${POT_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(pot SHARED ${POT_SRC})
ADD_LIBRARY(pot_stat STATIC ${POT_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: pot")
   target_precompile_headers(pot REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(pot_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(pot      cell linalg specialfunctions ${ext_libs})
TARGET_LINK_LIBRARIES(pot_stat cell_stat linalg_stat specialfunctions_stat ${ext_libs})


