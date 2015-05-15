#
#  Source files and headers in this directory
#
file(GLOB THERMOSTAT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB THERMOSTAT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${THERMOSTAT_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(thermostat SHARED ${THERMOSTAT_SRC})
ADD_LIBRARY(thermostat_stat STATIC ${THERMOSTAT_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(thermostat      mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(thermostat_stat mmath_stat ${ext_libs})


