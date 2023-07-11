#
#  Source files and headers in this directory
#
file(GLOB QOBJECTS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB QOBJECTS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${QOBJECTS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(qobjects SHARED ${QOBJECTS_SRC})
ADD_LIBRARY(qobjects_stat STATIC ${QOBJECTS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(qobjects      util molint linalg meigen specialfunctions ${ext_libs})
TARGET_LINK_LIBRARIES(qobjects_stat util_stat molint_stat linalg_stat meigen_stat specialfunctions_stat ${ext_libs})


