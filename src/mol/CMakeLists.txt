#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory atom...")
ADD_SUBDIRECTORY("atom")


#
#  Source files and headers in this directory
#
file(GLOB MOL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MOL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MOL_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(mol      SHARED ${MOL_SRC})
ADD_LIBRARY(mol_stat STATIC ${MOL_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(mol      atom_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(mol_stat atom_stat mmath_stat ${ext_libs})



