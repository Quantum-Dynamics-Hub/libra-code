#
#  Source files and headers in this directory
#
file(GLOB WFCGRID2_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB WFCGRID2_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${WFCGRID2_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(wfcgrid2 SHARED ${WFCGRID2_SRC})
ADD_LIBRARY(wfcgrid2_stat STATIC ${WFCGRID2_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in wfcgrid2")
   target_precompile_headers(wfcgrid2 REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(wfcgrid2_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(wfcgrid2      wfcgrid_stat nhamiltonian_stat linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(wfcgrid2_stat wfcgrid_stat nhamiltonian_stat linalg_stat random_stat meigen_stat specialfunctions_stat ${ext_libs})


