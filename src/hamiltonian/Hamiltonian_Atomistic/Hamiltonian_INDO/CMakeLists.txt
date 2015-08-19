#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_INDO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_INDO_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_INDO_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_indo SHARED ${HAMILTONIAN_INDO_SRC})
ADD_LIBRARY(hamiltonian_indo_stat STATIC ${HAMILTONIAN_INDO_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_indo      qchem_stat chemobjects_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES( hamiltonian_indo_stat qchem_stat chemobjects_stat mmath_stat ${ext_libs})


