#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_MODEL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_MODEL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_MODEL_HEADERS}) 


#set (EIGEN_INCLUDE ../../_Dependencies/Eigen/eigen-eigen-b9cd8366d4e8/ )
#INCLUDE_DIRECTORIES(${EIGEN_INCLUDE})



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_model SHARED ${HAMILTONIAN_MODEL_SRC})
ADD_LIBRARY(hamiltonian_model_stat STATIC ${HAMILTONIAN_MODEL_SRC})




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian_model       hamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(hamiltonian_model_stat  hamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})


