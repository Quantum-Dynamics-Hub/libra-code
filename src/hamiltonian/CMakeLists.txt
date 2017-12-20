#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_HEADERS}) 



MESSAGE("Going into subdirectory Hamiltonian_Generic...")
ADD_SUBDIRECTORY("Hamiltonian_Generic")

MESSAGE("Going into subdirectory nHamiltonian_Generic...")
ADD_SUBDIRECTORY("nHamiltonian_Generic")

MESSAGE("Going into subdirectory Hamiltonian_Model...")
ADD_SUBDIRECTORY("Hamiltonian_Model")


MESSAGE("Going into subdirectory Hamiltonian_Extern...")
ADD_SUBDIRECTORY("Hamiltonian_Extern")


MESSAGE("Going into subdirectory Hamiltonian_Atomistic...")
ADD_SUBDIRECTORY("Hamiltonian_Atomistic")






#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian SHARED ${HAMILTONIAN_SRC})
ADD_LIBRARY(hamiltonian_stat STATIC ${HAMILTONIAN_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian      
                      hamiltonian_atomistic_stat
                      hamiltonian_extern_stat 
                      hamiltonian_model_stat
                      hamiltonian_generic_stat 
                      nhamiltonian_generic_stat 
                      linalg_stat ${ext_libs}
                     )



TARGET_LINK_LIBRARIES(hamiltonian_stat
                      hamiltonian_atomistic_stat
                      hamiltonian_extern_stat 
                      hamiltonian_model_stat
                      hamiltonian_generic_stat
                      nhamiltonian_generic_stat
                      linalg_stat ${ext_libs}
                     )



