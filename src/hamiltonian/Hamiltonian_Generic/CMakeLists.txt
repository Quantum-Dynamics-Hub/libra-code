#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_GENERIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_GENERIC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_GENERIC_HEADERS}) 



INCLUDE_DIRECTORIES(/home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/)


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_generic SHARED ${HAMILTONIAN_GENERIC_SRC})
ADD_LIBRARY(hamiltonian_generic_stat STATIC ${HAMILTONIAN_GENERIC_SRC})




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian_generic      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(hamiltonian_generic_stat mmath_stat ${ext_libs})


