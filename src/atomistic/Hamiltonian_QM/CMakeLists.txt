#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_QM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_QM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_QM_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_qm SHARED ${HAMILTONIAN_QM_SRC})
ADD_LIBRARY(hamiltonian_qm_stat STATIC ${HAMILTONIAN_QM_SRC})



#if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
#   message(STATUS "Compiling using pre-compiled header support: in hamiltonian_qm")
#   target_precompile_headers(hamiltonian_qm REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
#   target_precompile_headers(hamiltonian_qm_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
#endif()




#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_qm    
                       basis_setups calculators
                       common_types
                       model_parameters control_parameters
                       qobjects chemobjects
                       linalg meigen specialfunctions)

TARGET_LINK_LIBRARIES( hamiltonian_qm_stat 
                       basis_setups_stat calculators_stat   
                       common_types_stat
                       model_parameters_stat control_parameters_stat 
                       qobjects_stat chemobjects_stat 
                       linalg_stat meigen_stat specialfunctions_stat)



