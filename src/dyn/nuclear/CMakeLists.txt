#
#  Source files and headers in this directory
#
file(GLOB NUCLEAR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB NUCLEAR_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${NUCLEAR_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(nuclear SHARED ${NUCLEAR_SRC})
ADD_LIBRARY(nuclear_stat STATIC ${NUCLEAR_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(nuclear      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(nuclear_stat linalg_stat ${ext_libs})


