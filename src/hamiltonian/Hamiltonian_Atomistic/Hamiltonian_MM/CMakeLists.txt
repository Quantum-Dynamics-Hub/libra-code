#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_MM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_MM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_MM_HEADERS}) 


MESSAGE("Going into subdirectory ForceField...")
ADD_SUBDIRECTORY("ForceField")


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_mm SHARED ${HAMILTONIAN_MM_SRC})
ADD_LIBRARY(hamiltonian_mm_stat STATIC ${HAMILTONIAN_MM_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_mm      forcefield_stat pot_stat cell_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES( hamiltonian_mm_stat forcefield_stat pot_stat cell_stat mmath_stat ${ext_libs})


