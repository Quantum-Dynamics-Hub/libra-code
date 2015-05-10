#
#  Source files and headers in this directory
#
file(GLOB OPERATORS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB OPERATORS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${OPERATORS_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(operators SHARED ${OPERATORS_SRC})
ADD_LIBRARY(operators_stat STATIC ${OPERATORS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(operators      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(operators_stat mmath_stat ${ext_libs})


