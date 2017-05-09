#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Source files and headers in this directory
#
file(GLOB QCHEM_TOOLS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB QCHEM_TOOLS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${QCHEM_TOOLS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(qchem_tools      SHARED ${QCHEM_TOOLS_SRC})
ADD_LIBRARY(qchem_tools_stat STATIC ${QCHEM_TOOLS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(qchem_tools      chemobjects_stat 
                      hamiltonian_stat qobjects_stat
                      basis_setups_stat control_parameters_stat
                      model_parameters_stat
                      linalg_stat ${ext_libs} )

TARGET_LINK_LIBRARIES(qchem_tools_stat      chemobjects_stat 
                      hamiltonian_stat qobjects_stat
                      basis_setups_stat control_parameters_stat
                      model_parameters_stat
                      linalg_stat ${ext_libs} )



