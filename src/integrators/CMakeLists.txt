#
#  Source files and headers in this directory
#
file(GLOB INTEGRATORS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB INTEGRATORS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${INTEGRATORS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(integrators      SHARED ${INTEGRATORS_SRC})
ADD_LIBRARY(integrators_stat STATIC ${INTEGRATORS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(integrators
                      linalg_stat meigen_stat util_stat specialfunctions_stat ${ext_libs} )

TARGET_LINK_LIBRARIES(integrators_stat
                      linalg_stat meigen_stat util_stat specialfunctions_stat ${ext_libs} )

