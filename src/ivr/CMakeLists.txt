#
#  Source files and headers in this directory
#
file(GLOB IVR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB IVR_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${IVR_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(ivr SHARED ${IVR_SRC})
ADD_LIBRARY(ivr_stat STATIC ${IVR_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: ivr")
   target_precompile_headers(ivr REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(ivr_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(ivr      meigen random linalg specialfunctions ${ext_libs})
TARGET_LINK_LIBRARIES(ivr_stat meigen_stat random_stat linalg_stat specialfunctions_stat ${ext_libs})


