#
#  Source files and headers in this directory
#
file(GLOB SYMMETRY_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SYMMETRY_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SYMMETRY_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(symmetry SHARED ${SYMMETRY_SRC})
ADD_LIBRARY(symmetry_stat STATIC ${SYMMETRY_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(symmetry      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(symmetry_stat linalg_stat ${ext_libs})


