#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_QM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_QM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_QM_HEADERS}) 



MESSAGE("Going into subdirectory Control_Parameters...")
ADD_SUBDIRECTORY("Control_Parameters")

MESSAGE("Going into subdirectory Model_Parameters...")
ADD_SUBDIRECTORY("Model_Parameters")



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_qm SHARED ${HAMILTONIAN_QM_SRC})
ADD_LIBRARY(hamiltonian_qm_stat STATIC ${HAMILTONIAN_QM_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_qm      calculators_stat model_parameters_stat control_parameters_stat qchem_stat chemobjects_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES( hamiltonian_qm_stat calculators_stat model_parameters_stat control_parameters_stat qchem_stat chemobjects_stat mmath_stat ${ext_libs})


