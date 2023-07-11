#
#  Source files and headers in this directory
#
file(GLOB SOLVERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SOLVERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SOLVERS_HEADERS}) 

# 
#SET( EIGEN_INCLUDE /home/Alexey_2/Soft/eigen-eigen-2249f9c22fe8/ )

set (EIGEN_INCLUDE ../math_meigen/ )
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE})

#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(solvers      SHARED ${SOLVERS_SRC})
ADD_LIBRARY(solvers_stat STATIC ${SOLVERS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(solvers      linalg ${ext_libs})
TARGET_LINK_LIBRARIES(solvers_stat linalg_stat ${ext_libs})

                                                 
