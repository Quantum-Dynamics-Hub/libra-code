#
#  Source files and headers in this directory
#
file(GLOB OPT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB OPT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${OPT_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(opt      SHARED ${OPT_SRC})
ADD_LIBRARY(opt_stat STATIC ${OPT_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(opt       linalg_stat  io_stat ${ext_libs} )
TARGET_LINK_LIBRARIES(opt_stat  linalg_stat  io_stat ${ext_libs} )




