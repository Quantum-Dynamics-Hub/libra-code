#
#  Source files and headers in this directory
#
file(GLOB FORCE_FIELD_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB FORCE_FIELD_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${FORCE_FIELD_HEADERS}) 



#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(forcefield SHARED ${FORCE_FIELD_SRC})
ADD_LIBRARY(forcefield_stat STATIC ${FORCE_FIELD_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( forcefield      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES( forcefield_stat mmath_stat ${ext_libs})


