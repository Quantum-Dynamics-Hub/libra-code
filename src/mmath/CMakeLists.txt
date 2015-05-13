MESSAGE("Going into subdirectory linalg...")
ADD_SUBDIRECTORY("linalg")

MESSAGE("Going into subdirectory specialfunctions...")
ADD_SUBDIRECTORY("specialfunctions")

MESSAGE("Going into subdirectory graph...")
ADD_SUBDIRECTORY("graph")

MESSAGE("Going into subdirectory operators...")
ADD_SUBDIRECTORY("operators")



#
#  Source files and headers in this directory
#
file(GLOB MMATH_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MMATH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MMATH_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(mmath      SHARED ${MMATH_SRC})
ADD_LIBRARY(mmath_stat STATIC ${MMATH_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(mmath      operators_stat graph_stat specialfunctions_stat linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(mmath_stat operators_stat graph_stat specialfunctions_stat linalg_stat ${ext_libs})

                                                 
