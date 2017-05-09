#
#  Source files and headers in this directory
#
file(GLOB UNIVERSE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB UNIVERSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${UNIVERSE_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(universe SHARED ${UNIVERSE_SRC})
ADD_LIBRARY(universe_stat STATIC ${UNIVERSE_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(universe      io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(universe_stat io_stat ${ext_libs})


