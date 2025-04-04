#
#  Source files and headers in this directory
#
file(GLOB sHAMILTONIAN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB sHAMILTONIAN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${sHAMILTONIAN_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(shamiltonian SHARED ${sHAMILTONIAN_SRC})
ADD_LIBRARY(shamiltonian_stat STATIC ${sHAMILTONIAN_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in shamiltonian")
   target_precompile_headers(shamiltonian REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(shamiltonian_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(shamiltonian_stat io_stat linalg_stat meigen_stat ${ext_libs})
TARGET_LINK_LIBRARIES(shamiltonian  shamiltonian_stat)


