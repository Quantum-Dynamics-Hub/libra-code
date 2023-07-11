#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_EXTERN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_EXTERN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_EXTERN_HEADERS}) 


#set (EIGEN_INCLUDE ../../_Dependencies/Eigen/eigen-eigen-b9cd8366d4e8/ )
#INCLUDE_DIRECTORIES(${EIGEN_INCLUDE})

#INCLUDE_DIRECTORIES(/home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/)


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_extern SHARED ${HAMILTONIAN_EXTERN_SRC})
ADD_LIBRARY(hamiltonian_extern_stat STATIC ${HAMILTONIAN_EXTERN_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in hamiltonian_extern")
   target_precompile_headers(hamiltonian_extern REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(hamiltonian_extern_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian_extern       hamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(hamiltonian_extern_stat  hamiltonian_generic_stat linalg_stat meigen_stat ${ext_libs})


