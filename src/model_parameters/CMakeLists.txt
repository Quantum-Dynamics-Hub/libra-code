#
#  Source files and headers in this directory
#
file(GLOB MODEL_PARAMETERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MODEL_PARAMETERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MODEL_PARAMETERS_HEADERS}) 



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(model_parameters SHARED ${MODEL_PARAMETERS_SRC})
ADD_LIBRARY(model_parameters_stat STATIC ${MODEL_PARAMETERS_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: model_parameters")
   target_precompile_headers(model_parameters REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(model_parameters_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( model_parameters      control_parameters qobjects   linalg ${ext_libs} )

TARGET_LINK_LIBRARIES( model_parameters_stat control_parameters_stat qobjects_stat linalg_stat ${ext_libs} )


