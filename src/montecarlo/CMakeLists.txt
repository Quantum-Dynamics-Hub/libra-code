#
#  Source files and headers in this directory
#
file(GLOB MC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MC_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(montecarlo      SHARED ${MC_SRC})
ADD_LIBRARY(montecarlo_stat STATIC ${MC_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(montecarlo      random_stat  linalg_stat  io_stat ${ext_libs} )
TARGET_LINK_LIBRARIES(montecarlo_stat random_stat  linalg_stat  io_stat ${ext_libs} )




