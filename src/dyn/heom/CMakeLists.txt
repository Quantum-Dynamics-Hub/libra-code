#
#  Source files and headers in this directory
#
file(GLOB HEOM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HEOM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HEOM_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(heom SHARED ${HEOM_SRC})
ADD_LIBRARY(heom_stat STATIC ${HEOM_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(heom      util_stat integrators_stat linalg_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(heom_stat util_stat integrators_stat linalg_stat specialfunctions_stat ${ext_libs})


