#
#  Source files and headers in this directory
#
file(GLOB nHAMILTONIAN_GENERIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB nHAMILTONIAN_GENERIC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${nHAMILTONIAN_GENERIC_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(nhamiltonian_generic SHARED ${nHAMILTONIAN_GENERIC_SRC})
ADD_LIBRARY(nhamiltonian_generic_stat STATIC ${nHAMILTONIAN_GENERIC_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(nhamiltonian_generic      linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(nhamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})


