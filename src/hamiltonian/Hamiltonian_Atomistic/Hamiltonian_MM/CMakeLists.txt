#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_MM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_MM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_MM_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_mm SHARED ${HAMILTONIAN_MM_SRC})
ADD_LIBRARY(hamiltonian_mm_stat STATIC ${HAMILTONIAN_MM_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_mm  
                       chemobjects_stat 
                       forcefield_stat
                       pot_stat
                       cell_stat 
                       linalg_stat specialfunctions_stat meigen_stat ${ext_libs})

TARGET_LINK_LIBRARIES( hamiltonian_mm_stat 
                       chemobjects_stat
                       forcefield_stat
                       pot_stat
                       cell_stat
                       linalg_stat specialfunctions_stat meigen_stat ${ext_libs})


