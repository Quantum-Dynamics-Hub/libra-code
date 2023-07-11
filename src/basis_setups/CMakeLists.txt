#
#  Source files and headers in this directory
#
file(GLOB BASIS_SETUPS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB BASIS_SETUPS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${BASIS_SETUPS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(basis_setups SHARED ${BASIS_SETUPS_SRC})
ADD_LIBRARY(basis_setups_stat STATIC ${BASIS_SETUPS_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: basis_setups")
   target_precompile_headers(basis_setups REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(basis_setups_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(basis_setups    basis  model_parameters qobjects linalg)

TARGET_LINK_LIBRARIES(basis_setups_stat basis_stat  model_parameters_stat qobjects_stat linalg_stat)


