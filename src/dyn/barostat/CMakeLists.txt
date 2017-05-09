#
#  Source files and headers in this directory
#
file(GLOB BAROSTAT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB BAROSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${BAROSTAT_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(barostat SHARED ${BAROSTAT_SRC})
ADD_LIBRARY(barostat_stat STATIC ${BAROSTAT_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(barostat      linalg_stat random_stat specialfunctions_stat io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(barostat_stat linalg_stat random_stat specialfunctions_stat io_stat ${ext_libs})


