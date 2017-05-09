#
#  Source files and headers in this directory
#
file(GLOB SPFUNC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SPFUNC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SPFUNC_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(specialfunctions      SHARED ${SPFUNC_SRC})
ADD_LIBRARY(specialfunctions_stat STATIC ${SPFUNC_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(specialfunctions      linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(specialfunctions_stat linalg_stat meigen_stat ${ext_libs})

                                                 
