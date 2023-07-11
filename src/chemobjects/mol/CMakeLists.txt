#
#  Source files and headers in this directory
#
file(GLOB MOL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MOL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MOL_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(mol SHARED ${MOL_SRC})
ADD_LIBRARY(mol_stat STATIC ${MOL_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: mol")
   target_precompile_headers(mol REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(mol_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(mol      universe rigidbody linalg)
TARGET_LINK_LIBRARIES(mol_stat universe_stat rigidbody_stat linalg_stat)


