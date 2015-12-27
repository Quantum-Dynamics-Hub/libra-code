#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_EXTERN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_EXTERN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_EXTERN_HEADERS}) 



INCLUDE_DIRECTORIES(/home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/)


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_extern SHARED ${HAMILTONIAN_EXTERN_SRC})
ADD_LIBRARY(hamiltonian_extern_stat STATIC ${HAMILTONIAN_EXTERN_SRC})




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian_extern       hamiltonian_generic_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(hamiltonian_extern_stat  hamiltonian_generic_stat mmath_stat ${ext_libs})


