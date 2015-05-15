#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_HEADERS}) 


#SET (p_eigen /home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/ )
#I=-I${p_eigen}


INCLUDE_DIRECTORIES(/home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/)


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian SHARED ${HAMILTONIAN_SRC})
ADD_LIBRARY(hamiltonian_stat STATIC ${HAMILTONIAN_SRC})




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(hamiltonian_stat mmath_stat ${ext_libs})


