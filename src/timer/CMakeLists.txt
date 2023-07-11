#
#  Source files and headers in this directory
#
file(GLOB TIMER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB TIMER_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${TIMER_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(timer SHARED ${TIMER_SRC})
ADD_LIBRARY(timer_stat STATIC ${TIMER_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(timer       ${ext_libs})
TARGET_LINK_LIBRARIES(timer_stat  ${ext_libs})


