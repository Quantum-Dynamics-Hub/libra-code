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



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( model_parameters      control_parameters_stat qobjects_stat
                       linalg_stat ${ext_libs} )

TARGET_LINK_LIBRARIES( model_parameters_stat control_parameters_stat qobjects_stat 
                       linalg_stat ${ext_libs} )


