#
#  Source files and headers in this directory
#
file(GLOB MOL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MOL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MOL_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(mol SHARED ${MOL_SRC})
ADD_LIBRARY(mol_stat STATIC ${MOL_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(mol      universe_stat rigidbody_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(mol_stat universe_stat rigidbody_stat mmath_stat ${ext_libs})


