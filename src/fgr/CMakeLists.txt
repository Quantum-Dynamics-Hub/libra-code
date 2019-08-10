#
#  Source files and headers in this directory
#
file(GLOB FGR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB FGR_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${FGR_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(fgr SHARED ${FGR_SRC})
ADD_LIBRARY(fgr_stat STATIC ${FGR_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(fgr      meigen_stat random_stat linalg_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(fgr_stat meigen_stat random_stat linalg_stat specialfunctions_stat ${ext_libs})


