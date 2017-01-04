#
#  Source files and headers in this directory
#
file(GLOB ANN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ANN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${ANN_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(ann_stat STATIC ${ANN_SRC})
ADD_LIBRARY(ann SHARED ${ANN_SRC})

#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(ann_stat specialfunctions_stat data_stat linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(ann      ann_stat)


