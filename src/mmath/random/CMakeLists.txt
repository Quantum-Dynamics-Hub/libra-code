#
#  Source files and headers in this directory
#
file(GLOB RANDOM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB RANDOM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${RANDOM_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(random SHARED ${RANDOM_SRC})
ADD_LIBRARY(random_stat STATIC ${RANDOM_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(random      ${ext_libs})
TARGET_LINK_LIBRARIES(random_stat ${ext_libs})


