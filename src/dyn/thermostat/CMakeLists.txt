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


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in thermostat")
   target_precompile_headers(thermostat REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(thermostat_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(thermostat      linalg_stat random_stat specialfunctions_stat io_stat ${ext_libs})
TARGET_LINK_LIBRARIES(thermostat_stat linalg_stat random_stat specialfunctions_stat io_stat ${ext_libs})


