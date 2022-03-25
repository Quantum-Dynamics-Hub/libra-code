#
#  Source files and headers in this directory
#
file(GLOB QTAG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB QTAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${QTAG_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(qtag SHARED ${QTAG_SRC})
ADD_LIBRARY(qtag_stat STATIC ${QTAG_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(qtag      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(qtag_stat linalg_stat ${ext_libs})


