#
#  Source files and headers in this directory
#
file(GLOB STATE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB STATE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${STATE_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(state SHARED ${STATE_SRC})
ADD_LIBRARY(state_stat STATIC ${STATE_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(state      atomistic_stat chemobjects_stat dyn_stat linalg_stat random_stat io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(state_stat atomistic_stat chemobjects_stat dyn_stat linalg_stat random_stat io_stat ${ext_libs})


