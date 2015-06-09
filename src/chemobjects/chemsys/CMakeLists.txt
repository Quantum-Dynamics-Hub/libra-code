#
#  Source files and headers in this directory
#
file(GLOB CHEMSYS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CHEMSYS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CHEMSYS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(chemsys SHARED ${CHEMSYS_SRC})
ADD_LIBRARY(chemsys_stat STATIC ${CHEMSYS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(chemsys      mol_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(chemsys_stat mol_stat mmath_stat ${ext_libs})


