#
#  Source files and headers in this directory
#
file(GLOB BASIS_SETUPS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB BASIS_SETUPS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${BASIS_SETUPS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(basis_setups SHARED ${BASIS_SETUPS_SRC})
ADD_LIBRARY(basis_setups_stat STATIC ${BASIS_SETUPS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(basis_setups      basis_stat qobjects_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(basis_setups_stat basis_stat qobjects_stat mmath_stat ${ext_libs})


