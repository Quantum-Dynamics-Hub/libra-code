#
#  Source files and headers in this directory
#
file(GLOB ELECTRONIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ELECTRONIC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${ELECTRONIC_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(electronic SHARED ${ELECTRONIC_SRC})
ADD_LIBRARY(electronic_stat STATIC ${ELECTRONIC_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(electronic      hamiltonian_stat 
                      linalg_stat meigen_stat random_stat operators_stat ${ext_libs} )

TARGET_LINK_LIBRARIES(electronic_stat hamiltonian_stat
                      linalg_stat meigen_stat random_stat operators_stat ${ext_libs} )


