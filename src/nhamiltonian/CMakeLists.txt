#
#  Source files and headers in this directory
#
file(GLOB nHAMILTONIAN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB nHAMILTONIAN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${nHAMILTONIAN_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(nhamiltonian SHARED ${nHAMILTONIAN_SRC})
ADD_LIBRARY(nhamiltonian_stat STATIC ${nHAMILTONIAN_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in nhamiltonian")
   target_precompile_headers(nhamiltonian REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(nhamiltonian_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(nhamiltonian_stat io_stat linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(nhamiltonian  nhamiltonian_stat)


