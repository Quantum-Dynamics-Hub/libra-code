#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_ATOMISTIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_ATOMISTIC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_ATOMISTIC_HEADERS}) 



#MESSAGE("Going into subdirectory Hamiltonian_Atomistic...")
#ADD_SUBDIRECTORY("Hamiltonian_Atomistic")


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_atomistic SHARED ${HAMILTONIAN_ATOMISTIC_SRC})
ADD_LIBRARY(hamiltonian_atomistic_stat STATIC ${HAMILTONIAN_ATOMISTIC_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_atomistic      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES( hamiltonian_atomistic_stat mmath_stat ${ext_libs})


