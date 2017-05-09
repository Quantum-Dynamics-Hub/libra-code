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



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_qm    
                       basis_setups_stat calculators_stat   
                       common_types_stat
                       model_parameters_stat control_parameters_stat 
                       qobjects_stat chemobjects_stat 
                       linalg_stat meigen_stat specialfunctions_stat ${ext_libs} )

TARGET_LINK_LIBRARIES( hamiltonian_qm_stat 
                       basis_setups_stat calculators_stat   
                       common_types_stat
                       model_parameters_stat control_parameters_stat 
                       qobjects_stat chemobjects_stat 
                       linalg_stat meigen_stat specialfunctions_stat ${ext_libs} )



