#
#  Source files and headers in this directory
#
file(GLOB SOLVERS_QM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SOLVERS_QM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SOLVERS_QM_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(solvers_qm      SHARED ${SOLVERS_QM_SRC})
ADD_LIBRARY(solvers_qm_stat STATIC ${SOLVERS_QM_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(solvers_qm      qobjects_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(solvers_qm_stat qobjects_stat mmath_stat ${ext_libs})

                                                 
