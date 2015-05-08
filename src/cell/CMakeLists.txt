#
#  Source files and headers in this directory
#
file(GLOB CELL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CELL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CELL_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(cell SHARED ${CELL_SRC})
ADD_LIBRARY(cell_stat STATIC ${CELL_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(cell      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(cell_stat mmath_stat ${ext_libs})


