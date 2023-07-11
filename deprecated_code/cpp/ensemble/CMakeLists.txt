#
#  Source files and headers in this directory
#
file(GLOB ENSEMBLE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ENSEMBLE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${ENSEMBLE_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(ensemble SHARED ${ENSEMBLE_SRC})
ADD_LIBRARY(ensemble_stat STATIC ${ENSEMBLE_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in ensemble")
   target_precompile_headers(ensemble REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(ensemble_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(ensemble      electronic_stat nuclear_stat hamiltonian_stat 
                      linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs} )

TARGET_LINK_LIBRARIES(ensemble_stat electronic_stat nuclear_stat hamiltonian_stat
                      linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs} )


