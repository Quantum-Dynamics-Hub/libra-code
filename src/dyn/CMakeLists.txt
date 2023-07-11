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

MESSAGE("Going into subdirectory wfcgrid2...")
ADD_SUBDIRECTORY("wfcgrid2")

#MESSAGE("Going into subdirectory ensemble...")
#ADD_SUBDIRECTORY("ensemble")

MESSAGE("Going into subdirectory gwp...")
ADD_SUBDIRECTORY("gwp")

MESSAGE("Going into subdirectory heom...")
ADD_SUBDIRECTORY("heom")

MESSAGE("Going into subdirectory qtag...")
ADD_SUBDIRECTORY("qtag")


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



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in dyn")
   target_precompile_headers(dyn REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(dyn_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(dyn      
                      heom
                      qtag
                      gwp
                      wfcgrid
                      wfcgrid2
                      barostat thermostat
                      electronic 
                      nuclear
                      nhamiltonian
                      converters_stat
                      linalg meigen util specialfunctions ${ext_libs} )

TARGET_LINK_LIBRARIES(dyn_stat
                      heom_stat
                      qtag_stat
                      gwp_stat
                      wfcgrid_stat 
                      wfcgrid2_stat 
                      barostat_stat thermostat_stat 
                      electronic_stat 
                      nuclear_stat 
                      nhamiltonian_stat
                      converters_stat
                      linalg_stat meigen_stat util_stat specialfunctions_stat ${ext_libs} )




