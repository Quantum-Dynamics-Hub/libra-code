#
#  Source files and headers in this directory
#
file(GLOB LINALG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB LINALG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${LINALG_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(linalg      SHARED ${LINALG_SRC})
ADD_LIBRARY(linalg_stat STATIC ${LINALG_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(linalg      io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(linalg_stat io_stat ${ext_libs})

                                                 
