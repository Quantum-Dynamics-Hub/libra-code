#
#  Source files and headers in this directory
#
file(GLOB UTIL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB UTIL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${UTIL_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(util      SHARED ${UTIL_SRC})
ADD_LIBRARY(util_stat STATIC ${UTIL_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(util      ${ext_libs})
TARGET_LINK_LIBRARIES(util_stat ${ext_libs})

                                                 
