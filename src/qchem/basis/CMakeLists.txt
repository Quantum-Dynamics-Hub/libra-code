#
#  Source files and headers in this directory
#
file(GLOB BASIS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB BASIS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${BASIS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(basis SHARED ${BASIS_SRC})
ADD_LIBRARY(basis_stat STATIC ${BASIS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(basis      qobjects_stat molint_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(basis_stat qobjects_stat molint_stat mmath_stat ${ext_libs})


