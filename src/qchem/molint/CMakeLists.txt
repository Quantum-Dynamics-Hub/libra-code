#
#  Source files and headers in this directory
#
file(GLOB MOLINT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MOLINT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MOLINT_HEADERS}) 

#set ( MOLINT_HEADERS libmolint.h Overlaps.h )
#set ( MOLINT_SRC libmolint.cpp Overlaps.cpp ${MOLINT} )


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(molint SHARED ${MOLINT_SRC})
ADD_LIBRARY(molint_stat STATIC ${MOLINT_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(molint      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(molint_stat mmath_stat ${ext_libs})


