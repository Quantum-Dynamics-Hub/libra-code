#
#  Source files and headers in this directory
#
file(GLOB CONTEXT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CONTEXT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CONTEXT_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(context      SHARED ${CONTEXT_SRC})
ADD_LIBRARY(context_stat STATIC ${CONTEXT_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(context      linalg_stat io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(context_stat linalg_stat io_stat ${ext_libs})

                                                 
