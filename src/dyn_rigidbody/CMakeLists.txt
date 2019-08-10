#
#  Source files and headers in this directory
#
file(GLOB RB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB RB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${RB_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(rigidbody SHARED ${RB_SRC})
ADD_LIBRARY(rigidbody_stat STATIC ${RB_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(rigidbody      linalg_stat meigen_stat specialfunctions_stat io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(rigidbody_stat linalg_stat meigen_stat specialfunctions_stat io_stat ${ext_libs})


