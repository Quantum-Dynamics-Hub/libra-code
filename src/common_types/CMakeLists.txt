#
#  Source files and headers in this directory
#
file(GLOB COMMON_TYPES_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB COMMON_TYPES_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${COMMON_TYPES_HEADERS}) 



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(common_types  SHARED ${COMMON_TYPES_SRC})
ADD_LIBRARY(common_types_stat STATIC ${COMMON_TYPES_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( common_types       ${ext_libs})
TARGET_LINK_LIBRARIES( common_types_stat  ${ext_libs})


