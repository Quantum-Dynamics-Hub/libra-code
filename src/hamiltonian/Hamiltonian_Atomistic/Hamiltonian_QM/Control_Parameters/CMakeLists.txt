#
#  Source files and headers in this directory
#
file(GLOB CONTROL_PARAMETERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CONTROL_PARAMETERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CONTROL_PARAMETERS_HEADERS}) 



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(control_parameters SHARED ${CONTROL_PARAMETERS_SRC})
ADD_LIBRARY(control_parameters_stat STATIC ${CONTROL_PARAMETERS_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( control_parameters      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES( control_parameters_stat mmath_stat ${ext_libs})


