#
#  Source files and headers in this directory
#
file(GLOB GWP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB GWP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${GWP_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(gwp SHARED ${GWP_SRC})
ADD_LIBRARY(gwp_stat STATIC ${GWP_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(gwp      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(gwp_stat mmath_stat ${ext_libs})


