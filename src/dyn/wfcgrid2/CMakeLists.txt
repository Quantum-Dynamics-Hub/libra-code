#
#  Source files and headers in this directory
#
file(GLOB WFCGRID2_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB WFCGRID2_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${WFCGRID2_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(wfcgrid2 SHARED ${WFCGRID2_SRC})
ADD_LIBRARY(wfcgrid2_stat STATIC ${WFCGRID2_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(wfcgrid2      wfcgrid_stat hamiltonian_stat linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(wfcgrid2_stat wfcgrid_stat hamiltonian_stat linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs})


