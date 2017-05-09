#
#  Source files and headers in this directory
#
file(GLOB POT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB POT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${POT_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(pot SHARED ${POT_SRC})
ADD_LIBRARY(pot_stat STATIC ${POT_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(pot      cell_stat linalg_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(pot_stat cell_stat linalg_stat specialfunctions_stat ${ext_libs})


