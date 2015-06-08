MESSAGE("Going into subdirectory linalg...")
ADD_SUBDIRECTORY("linalg")

MESSAGE("Going into subdirectory specialfunctions...")
ADD_SUBDIRECTORY("specialfunctions")

MESSAGE("Going into subdirectory graph...")
ADD_SUBDIRECTORY("graph")

MESSAGE("Going into subdirectory operators...")
ADD_SUBDIRECTORY("operators")

MESSAGE("Going into subdirectory random...")
ADD_SUBDIRECTORY("random")

MESSAGE("Going into subdirectory data...")
ADD_SUBDIRECTORY("data")

MESSAGE("Going into subdirectory ann...")
ADD_SUBDIRECTORY("ann")

MESSAGE("Going into subdirectory meigen...")
ADD_SUBDIRECTORY("meigen")

MESSAGE("Going into subdirectory symmetry...")
ADD_SUBDIRECTORY("symmetry")



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
TARGET_LINK_LIBRARIES(mmath      symmetry_stat meigen_stat ann_stat data_stat random_stat operators_stat graph_stat specialfunctions_stat linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(mmath_stat symmetry_stat meigen_stat ann_stat data_stat random_stat operators_stat graph_stat specialfunctions_stat linalg_stat ${ext_libs})

                                                 
