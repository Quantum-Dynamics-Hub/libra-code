#
#  Source files and headers in this directory
#
file(GLOB DATA_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB DATA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${DATA_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(data SHARED ${DATA_SRC})
ADD_LIBRARY(data_stat STATIC ${DATA_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(data      ${ext_libs})
TARGET_LINK_LIBRARIES(data_stat ${ext_libs})


