#
#  Source files and headers in this directory
#
file(GLOB WFCGRID_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB WFCGRID_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${WFCGRID_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(wfcgrid SHARED ${WFCGRID_SRC})
ADD_LIBRARY(wfcgrid_stat STATIC ${WFCGRID_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(wfcgrid      hamiltonian_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(wfcgrid_stat hamiltonian_stat mmath_stat ${ext_libs})


