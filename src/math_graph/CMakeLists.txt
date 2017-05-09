#
#  Source files and headers in this directory
#
file(GLOB GRAPH_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB GRAPH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${GRAPH_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(graph      SHARED ${GRAPH_SRC})
ADD_LIBRARY(graph_stat STATIC ${GRAPH_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(graph      ${ext_libs})
TARGET_LINK_LIBRARIES(graph_stat ${ext_libs})

                                                 
