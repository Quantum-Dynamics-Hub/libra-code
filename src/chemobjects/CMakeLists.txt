#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory universe...")
ADD_SUBDIRECTORY("universe")

MESSAGE("Going into subdirectory mol...")
ADD_SUBDIRECTORY("mol")

MESSAGE("Going into subdirectory chemsys...")
ADD_SUBDIRECTORY("chemsys")




#
#  Source files and headers in this directory
#
file(GLOB CHEMOBJECTS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CHEMOBJECTS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CHEMOBJECTS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(chemobjects      SHARED ${CHEMOBJECTS_SRC})
ADD_LIBRARY(chemobjects_stat STATIC ${CHEMOBJECTS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(chemobjects      chemsys_stat mol_stat universe_stat linalg_stat ${ext_libs})   
TARGET_LINK_LIBRARIES(chemobjects_stat chemsys_stat mol_stat universe_stat linalg_stat ${ext_libs})  



