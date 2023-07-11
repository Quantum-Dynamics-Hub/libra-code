#
#  Source files and headers in this directory
#
file(GLOB FGR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB FGR_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${FGR_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(fgr SHARED ${FGR_SRC})
ADD_LIBRARY(fgr_stat STATIC ${FGR_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: fgr")
   target_precompile_headers(fgr REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(fgr_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(fgr      meigen random linalg specialfunctions)
TARGET_LINK_LIBRARIES(fgr_stat meigen_stat random_stat linalg_stat specialfunctions_stat)


