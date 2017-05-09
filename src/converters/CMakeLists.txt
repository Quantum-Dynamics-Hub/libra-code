#
#  Source files and headers in this directory
#
file(GLOB CONVERTERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CONVERTERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CONVERTERS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(converters      SHARED ${CONVERTERS_SRC})
ADD_LIBRARY(converters_stat STATIC ${CONVERTERS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(converters      dyn_stat chemobjects_stat linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(converters_stat dyn_stat chemobjects_stat linalg_stat ${ext_libs})



