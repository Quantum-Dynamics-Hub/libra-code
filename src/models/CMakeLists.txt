
#
#  Source files and headers in this directory
#
file(GLOB MODELS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MODELS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MODELS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(models SHARED ${MODELS_SRC})
ADD_LIBRARY(models_stat STATIC ${MODELS_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in models")
   target_precompile_headers(models REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(models_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(models       linalg meigen ${ext_libs})
TARGET_LINK_LIBRARIES(models_stat  linalg_stat meigen_stat ${ext_libs})


