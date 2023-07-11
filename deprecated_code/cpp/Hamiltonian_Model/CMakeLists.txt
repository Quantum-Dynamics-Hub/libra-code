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



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in hamiltonian_model")
   target_precompile_headers(hamiltonian_model REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(hamiltonian_model_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian_model       hamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(hamiltonian_model_stat  hamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})


