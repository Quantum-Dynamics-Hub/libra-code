#
#  Source files and headers in this directory
#
file(GLOB ENSEMBLE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ENSEMBLE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${ENSEMBLE_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(ensemble SHARED ${ENSEMBLE_SRC})
ADD_LIBRARY(ensemble_stat STATIC ${ENSEMBLE_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(ensemble      electronic_stat nuclear_stat hamiltonian_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(ensemble_stat electronic_stat nuclear_stat hamiltonian_stat mmath_stat ${ext_libs})


