#
#  Building sub-directories
#
MESSAGE("Going into subdirectory nuclear...")
ADD_SUBDIRECTORY("nuclear")

#MESSAGE("Going into subdirectory rigidbody...")
#ADD_SUBDIRECTORY("rigidbody")

MESSAGE("Going into subdirectory electronic...")
ADD_SUBDIRECTORY("electronic")

MESSAGE("Going into subdirectory thermostat...")
ADD_SUBDIRECTORY("thermostat")

MESSAGE("Going into subdirectory barostat...")
ADD_SUBDIRECTORY("barostat")

MESSAGE("Going into subdirectory wfcgrid...")
ADD_SUBDIRECTORY("wfcgrid")

MESSAGE("Going into subdirectory ensemble...")
ADD_SUBDIRECTORY("ensemble")

MESSAGE("Going into subdirectory gwp...")
ADD_SUBDIRECTORY("gwp")


#
#  Source files and headers in this directory
#
file(GLOB DYN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB DYN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${DYN_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(dyn      SHARED ${DYN_SRC})
ADD_LIBRARY(dyn_stat STATIC ${DYN_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(dyn      gwp_stat ensemble_stat 
                      wfcgrid_stat 
                      barostat_stat thermostat_stat 
                      electronic_stat 
                      nuclear_stat 
                      hamiltonian_stat
                      linalg_stat meigen_stat specialfunctions_stat ${ext_libs} )

TARGET_LINK_LIBRARIES(dyn_stat      gwp_stat ensemble_stat 
                      wfcgrid_stat 
                      barostat_stat thermostat_stat 
                      electronic_stat 
                      nuclear_stat 
                      hamiltonian_stat
                      linalg_stat meigen_stat specialfunctions_stat ${ext_libs} )




